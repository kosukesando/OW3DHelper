module OW3DHelper

# using SpecialFunctions
using Reexport
using Printf
using FFTW
using OhMyThreads: @tasks, @local, @set
using NCDatasets
using Statistics: mean
using Printf
using DataStructures: OrderedDict
using DSP
using Serialization
using MLStyle
using CairoMakie
using SparseArrays


include("./Structs.jl")
include("./Constants.jl")
include("./IO.jl")
include("./Post-processing.jl")
include("./Video.jl")

export
    OW3DInput,
    EPFile,
    JSpec,
    KinematicSetting,
    calc_etaphi,
    # calc_etaphi_sec,
    export_ow3d_inp,
    export_ow3d_init,
    export_nc_ep,
    export_nc_hilbert,
    export_nc_4phase,
    plot_raw2d,
    plot_raw3d,
    plot_4p,
    getEPlist,
    open_EP,
    open_Kinematics,
    generate_init

function mcallister_mwd(foverfp::Float64, skewd::Float64)
    # Based on parameterisation in McAllister(2019)
    if foverfp <= 1
        theta = 0.0
    else
        theta = (tanh(2 * foverfp - 3.5) / tanh(1.5) + 1) * skewd / 2
    end
    return theta
end

function calc_k_omega_a(oi)
    kminx = 2 * pi / (oi.dx * oi.nx)
    kminy = 2 * pi / (oi.dy * oi.ny)

    kxi = kminx * (-(oi.nx - 1)/2:1:(oi.nx-1)/2)
    kyi = kminy * (-(oi.ny - 1)/2:1:(oi.ny-1)/2)
    # max of kxi converges to pi/oi.dx with nx->inf quite quickly
    kx = kxi[abs.(kxi).<=oi.kmaxx] # wavenumber components in x-direction
    ky = kyi[abs.(kyi).<=oi.kmaxy] # wavenumber components in y-direction

    #Note: wave components calculated directly from the power spectrum as a
    #characterisitc feature of the NewWave formulation but it would usually be
    #necessary to transform the power spectrum into an amplitude spectrum
    #before calculating the wave components
    nkx = length(kx)
    nky = length(ky)
    kmatg = zeros(nkx, nky)
    ωmatg = zeros(nkx, nky)
    fmatg = zeros(nkx, nky)
    ampg = zeros(nkx, nky)
    dirg = zeros(nkx, nky)

    # Generate spectrum. Type J is JONSWAP spectrum
    fm = 0.0836708928777909 #depth = 200m
    if typeof(oi.spec) == JSpec

        #NOTE: the value of fm is depth dependent!
        #peak of JONSWAP ensures kp of 0.02796 [m^-1]
        #fm=0.0787103953986135; #depth = 50m

        for i = 1:nkx
            for j = 1:nky
                kAbs = sqrt((kx[i])^2 + (ky[j])^2) #wavenumber of component
                theta = atan(ky[j], kx[i]) * 180 / pi #direction of component
                ω = sqrt(kAbs * g .* tanh(kAbs * oi.depth)) #natural frequency of component
                f = (1 / (2 * pi)) * ω #natural frequencies in Hertz
                kmatg[i, j] = kAbs
                dirg[i, j] = theta
                ωmatg[i, j] = ω
                fmatg[i, j] = f
            end
        end

        sj = zeros(nkx, nky) #bandwidth parameters of JONSWAP
        sj[fmatg.<=fm] .= 0.07
        sj[fmatg.>fm] .= 0.09

        #JONSWAP spectrum based on wavenumber
        kd = oi.k0 * oi.depth
        D = (g^2) * ((2 * pi)^(-6)) * (0.5 * g) * (tanh(kd) + kd * ((sech(kd))^2))
        h1 = fmatg .^ (-6) #differs from frequency version due to Jacobian
        h2 = exp.((-5 / 4) * ((fmatg / fm) .^ (-4))) #exponential function
        h3 = oi.spec.γ .^ (exp.(-0.5 * ((((fmatg / fm) .- 1) ./ (sj)) .^ 2))) #peak enhancement
        Spec = D * h1 .* h2 .* h3
        # Ignore "alpha" coefficient since this is a scaling parameter
        # Distribution valid for S(ω) & S(f) since scaled with normalisation

        if oi.spreading_type == "W"
            Sdir = exp.(-0.5 * (dirg .^ 2) ./ (oi.spreading_param^2)) # A wrapped normal spreading function with parameter sigma
        end
        ampg_newwave = 1 ./ kmatg .* Spec .* Sdir #amplitude of component - S(kx,ky)
        replace!(ampg_newwave, NaN => 0)
    else
        println("spectrum not recognised")
    end

    if oi.twist_type == "mcallister"
        dirg = dirg + mcallister_mwd.(fmatg / fm, oi.twist_angle) # MWD as a function of freq
    end
    dirg = dirg .+ oi.mwd

    kstarth = sqrt(oi.kmaxx * oi.kmaxy) - 1.00 * oi.k0
    kmin = sqrt(kminx * kminy)
    kmax = sqrt(oi.kmaxx * oi.kmaxy)
    kstartl = 0.5 * oi.k0
    Ism = abs.(kmatg) .>= kstarth .&& abs.(kmatg) .<= kmax
    Isml = abs.(kmatg) .>= kmin .&& abs.(kmatg) .<= kstartl
    ampg_newwave[Ism] .= ampg_newwave[Ism] .* (0.5 * (cos.(((abs.(kmatg[Ism]) .- kstarth) ./ (kmax .- kstarth)) * pi) .+ 1))
    ampg_newwave[Isml] .= ampg_newwave[Isml] .* (0.5 * (cos.(((abs.(kmatg[Isml]) .- kstartl) ./ (kmin .- kstartl)) * pi) .+ 1))
    ampg_newwave[abs.(kmatg).>=kmax] .= 0
    ampg_newwave[abs.(kmatg).<=kmin] .= 0
    ampg_newwave[isnan.(ampg_newwave)] .= 0

    kxmatg = kmatg .* cos.(dirg * (pi / 180)) #Component of k in x direction
    kymatg = kmatg .* sin.(dirg * (pi / 180)) #Component of k in y direction

    ampg_newwave_norm = oi.A * (ampg_newwave ./ sum(ampg_newwave)) #  normalise to control size of event

    kxmatg, kymatg, ωmatg, ampg_newwave_norm, dirg
end

function calc_etaphi(oi::OW3DInput, t::Number; order=1)
    kxmatg, kymatg, ωmatg, ampg_newwave_norm, dirg = calc_k_omega_a(oi)
    η = calc_eta(oi, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm)
    ϕ = calc_phi(oi, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm, η)
    if order == 1
        return η, ϕ
    elseif order == 2
        η₂₂, ϕ₂₂, η₂₀, ϕ₂₀ = calc_etaphi_sec(oi, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm, dirg)
        return η, ϕ, η₂₂, ϕ₂₂, η₂₀, ϕ₂₀
    end
end

function calc_etaphi(oi::OW3DInput, ts::Int, te::Int, x, y)
    kxmatg, kymatg, ωmatg, ampg_newwave_norm, dirg = calc_k_omega_a(oi)
    t = ts:te
    η = calc_eta(oi, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm, x, y)
    # ϕ = calc_phi(oi, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm, x, y)
    ϕ = zeros(length(t))
    η, ϕ
end

function calc_eta(oi, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm)
    (nkx, nky) = size(ampg_newwave_norm)
    # Generate full-domain
    xvec = oi.dx * (-(oi.nx - 1)/2:1:(oi.nx-1)/2)
    yvec = oi.dy * (-(oi.ny - 1)/2:1:(oi.ny-1)/2)
    xmat = xvec .* ones(length(yvec))'
    ymat = yvec' .* ones(length(xvec))
    # Calculate linear free surface
    println("Calculating linear free surface(omtt)")
    η = @tasks for kij = eachindex(1:nky*nkx)
        @set reducer = .+
        local ki = 1 + (kij - 1) % nky
        local kj = 1 + (kij - 1) ÷ nky
        @local η_kj = zeros(oi.nx, oi.ny)
        local phasei = kxmatg[ki, kj] .* xmat .+ kymatg[ki, kj] .* ymat .- ωmatg[ki, kj] * t .+ deg2rad(oi.ϕ)
        ampg_newwave_norm[ki, kj] .* cos.(phasei)
    end
    return η
end

function calc_eta(oi, kxmatg, kymatg, ωmatg, t_vec, ampg_newwave_norm, x::Float64, y::Float64)
    (nkx, nky) = size(ampg_newwave_norm)
    η = zeros(length(t_vec))
    for (i, t) in enumerate(t_vec)
        for kj = eachindex(1:nky), ki = eachindex(1:nkx)
            kx = kxmatg[ki, kj]
            ky = kymatg[ki, kj]
            ω = ωmatg[ki, kj]
            an = ampg_newwave_norm[ki, kj]
            phasei = kx * x + ky * y - ω * t + deg2rad(oi.ϕ)
            etacomp = an * cos(phasei)
            η[i] += etacomp
        end
    end
    η
end

function calc_phi(oi, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm, η)
    (nkx, nky) = size(ampg_newwave_norm)
    # Generate full-domain
    xvec = oi.dx * (-(oi.nx - 1)/2:1:(oi.nx-1)/2)
    yvec = oi.dy * (-(oi.ny - 1)/2:1:(oi.ny-1)/2)
    xmat = xvec .* ones(length(yvec))'
    ymat = yvec' .* ones(length(xvec))
    kmatg = sqrt.(kxmatg .^ 2 .+ kymatg .^ 2)
    println("Calculating velocity potential at free surface")
    ϕ = @tasks for kij = eachindex(1:nky*nkx)
        @set reducer = .+
        local ki = 1 + (kij - 1) % nky
        local kj = 1 + (kij - 1) ÷ nky
        @local ϕ_kj = zeros(oi.nx, oi.ny)
        local phasei = kxmatg[ki, kj] .* xmat .+ kymatg[ki, kj] .* ymat .- ωmatg[ki, kj] * t .+ deg2rad(oi.ϕ)
        (ampg_newwave_norm[ki, kj] ./ (ωmatg[ki, kj] .+ 0.000000001)) .* g .* ((cosh.(kmatg[ki, kj] .* (η .+ oi.depth)) ./ cosh.(kmatg[ki, kj] .* oi.depth)) .* sin.(phasei))
    end
    ϕ
end

function calc_etaphi_origin(oi::OW3DInput, ts::Int, te::Int)
    calc_etaphi(oi, ts, te, 0.0, 0.0)
end

function calc_etaphi_sec(oi::OW3DInput, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm, dirg)
    nx = oi.nx
    ny = oi.ny
    phase = deg2rad(oi.ϕ)
    kminx = 2 * pi / (oi.dx * nx)
    kminy = 2 * pi / (oi.dy * ny)
    kxi = kminx * (-(nx - 1)/2:1:(nx-1)/2)
    kyi = kminy * (-(ny - 1)/2:1:(ny-1)/2)
    # max of kxi converges to pi/oi.dx with nx->inf quite quickly
    kx = kxi[abs.(kxi).<=oi.kmaxx] # wavenumber components in x-direction
    ky = kyi[abs.(kyi).<=oi.kmaxy] # wavenumber components in y-direction

    #Note: wave components calculated directly from the power spectrum as a
    #characterisitc feature of the NewWave formulation but it would usually be
    #necessary to transform the power spectrum into an amplitude spectrum
    #before calculating the wave components
    nkx = length(kx)
    nky = length(ky)
    depth = oi.depth
    kmatg = sqrt.(kxmatg .^ 2 .+ kymatg .^ 2)

    C22 = zeros(ComplexF64, 2 * (nx - 1) + 1, 2 * (ny - 1) + 1)
    V22 = zeros(ComplexF64, 2 * (nx - 1) + 1, 2 * (ny - 1) + 1)
    C20 = zeros(ComplexF64, 2 * (nx - 1) + 1, 2 * (ny - 1) + 1)
    V20 = zeros(ComplexF64, 2 * (nx - 1) + 1, 2 * (ny - 1) + 1)

    ath = 10^(-6) # amplitude threshold for calculation of second order interactions

    #First calculate self interaction term.

    println("Calculating self interaction term")
    maxes = 0

    # for q00=1:1:length(iia)
    for i = 1:nkx
        for j = 1:nky
            if ampg_newwave_norm[i, j] < ath
                continue
            end
            e2sx = minimum(abs.(kxi .- (kxmatg[i, j] + kxmatg[i, j])))
            e2sy = minimum(abs.(kyi .- (kymatg[i, j] + kymatg[i, j])))
            iix2s = argmin(abs.(kxi .- (kxmatg[i, j] + kxmatg[i, j])))
            iiy2s = argmin(abs.(kyi .- (kymatg[i, j] + kymatg[i, j])))

            maxes = maximum([maxes e2sx e2sy])

            if maxes > (10^-12)
                "warning: misassigned component"
            end

            as2s = ((ampg_newwave_norm[i, j]^2) * kmatg[i, j] / (4 * tanh(kmatg[i, j] * depth))) * (2 + 3 / ((sinh(kmatg[i, j] * depth))^2))
            av2s = ((ampg_newwave_norm[i, j]^2) * (3 * ωmatg[i, j] / 8) * cosh(2 * kmatg[i, j] * depth) / ((sinh(kmatg[i, j] * depth))^4))
            p2s = (phase - t * ωmatg[i, j]) + (phase - t * ωmatg[i, j])
            c2s = complex(as2s * cos(p2s), as2s * sin(p2s))
            v2s = complex(av2s * cos(p2s - pi / 2), av2s * sin(p2s - pi / 2))

            C22[iix2s, iiy2s] = C22[iix2s, iiy2s] + c2s
            V22[iix2s, iiy2s] = V22[iix2s, iiy2s] + v2s
        end
    end

    #Now calculate the interaction between each component

    println("Calculating second order cross-interactions...")

    maxec = 0

    for ij1 = 1:nkx*nky-1
        i1 = 1 + (ij1 - 1) % nky
        j1 = 1 + (ij1 - 1) ÷ nky
        sinh1squared = (sinh(kmatg[i1, j1] * depth))^2
        tanh1 = tanh(kmatg[i1, j1] * depth)
        for ij2 = ij1+1:nkx*nky
            i2 = 1 + (ij2 - 1) % nky
            j2 = 1 + (ij2 - 1) ÷ nky
            sinh2squared = (sinh(kmatg[i2, j2] * depth))^2
            tanh2 = tanh(kmatg[i2, j2] * depth)
            cosdt = cosd(dirg[i1, j1] - dirg[i2, j2]) # cos of the angle between the constituents
            pmodk = sqrt((kxmatg[i1, j1] + kxmatg[i2, j2])^2 + (kymatg[i1, j1] + kymatg[i2, j2])^2)  #mod(k1+k2)
            mmodk = sqrt((kxmatg[i1, j1] - kxmatg[i2, j2])^2 + (kymatg[i1, j1] - kymatg[i2, j2])^2)  #mod(k1-k2)
            Dp = (ωmatg[i1, j1] + ωmatg[i2, j2])^2 - g * pmodk * tanh(pmodk * depth) #Equation 21 in Dalzell
            Dm = (ωmatg[i1, j1] - ωmatg[i2, j2])^2 - g * mmodk * tanh(mmodk * depth) #Equation 22 in Dalzell
            tanhs = tanh1 * tanh2
            if tanhs < 0
                println("Alert! tanhs=$(tanhs)")
            end
            Ap1 = -(ωmatg[i1, j1] * ωmatg[i2, j2] * (ωmatg[i1, j1] + ωmatg[i2, j2]) / Dp) * (1 - cosdt / tanhs) #First term in eq 17 of Dalzell
            Ap2 = (1 / (2 * Dp)) * ((ωmatg[i1, j1]^3) / sinh1squared + (ωmatg[i2, j2]^3) / sinh2squared) #Second term in eq 17 of Dalzell
            Ap = Ap1 + Ap2
            Am1 = (ωmatg[i1, j1] * ωmatg[i2, j2] * (ωmatg[i1, j1] - ωmatg[i2, j2]) / Dm) * (1 + cosdt / tanhs) #First term in eq 18 of Dalzell
            Am2 = (1 / (2 * Dm)) * ((ωmatg[i1, j1]^3) / sinh1squared - (ωmatg[i2, j2]^3) / sinh2squared) #Second term in eq 18 of Dalzell
            Am = Am1 + Am2
            Bp0 = (ωmatg[i1, j1]^2 + ωmatg[i2, j2]^2) / (2 * g) #First term in 19 and 20 of Dalzell
            Bp1 = -((ωmatg[i1, j1] * ωmatg[i2, j2]) / (2 * g)) * (1 - cosdt / tanhs) * ((ωmatg[i1, j1] + ωmatg[i2, j2])^2 + g * pmodk * tanh(pmodk * depth)) / Dp #Second term in eq 19
            Bp2 = ((ωmatg[i1, j1] + ωmatg[i2, j2]) / (2 * g * Dp)) * (((ωmatg[i1, j1]^3) / sinh1squared) + ((ωmatg[i2, j2]^3) / sinh2squared)) #Third term in eq 19
            Bm1 = ((ωmatg[i1, j1] * ωmatg[i2, j2]) / (2 * g)) * (1 + cosdt / tanhs) * ((ωmatg[i1, j1] - ωmatg[i2, j2])^2 + g * mmodk * tanh(mmodk * depth)) / Dm #Second term in eq 20
            Bm2 = ((ωmatg[i1, j1] - ωmatg[i2, j2]) / (2 * g * Dm)) * (((ωmatg[i1, j1]^3) / sinh1squared) - ((ωmatg[i2, j2]^3) / sinh2squared)) #Third term in eq 19
            Bp = Bp0 + Bp1 + Bp2
            Bm = Bp0 + Bm1 + Bm2

            e2px = minimum(abs.(kxi .- (kxmatg[i1, j1] + kxmatg[i2, j2])))
            e2py = minimum(abs.(kyi .- (kymatg[i1, j1] + kymatg[i2, j2])))
            e2mx = minimum(abs.(kxi .- (kxmatg[i2, j2] - kxmatg[i1, j1])))
            e2my = minimum(abs.(kyi .- (kymatg[i2, j2] - kymatg[i1, j1])))
            iix2p = argmin(abs.(kxi .- (kxmatg[i1, j1] + kxmatg[i2, j2])))
            iiy2p = argmin(abs.(kyi .- (kymatg[i1, j1] + kymatg[i2, j2])))
            iix2m = argmin(abs.(kxi .- (kxmatg[i2, j2] - kxmatg[i1, j1])))
            iiy2m = argmin(abs.(kyi .- (kymatg[i2, j2] - kymatg[i1, j1])))

            maxec = maximum([maxec, e2px, e2py, e2mx, e2my])

            if maxec > (10^-12)
                println("warning: misassigned component")
            end

            a22 = ampg_newwave_norm[i1, j1] * ampg_newwave_norm[i2, j2]
            p22 = (phase - t * ωmatg[i1, j1]) + (phase - t * ωmatg[i2, j2])
            c22 = complex(a22 * Bp * cos(p22), a22 * Bp * sin(p22))
            v22 = complex(a22 * Ap * cos(p22 - pi / 2), a22 * Ap * sin(p22 - pi / 2))

            a20 = ampg_newwave_norm[i1, j1] * ampg_newwave_norm[i2, j2]
            p20 = (phase - t * ωmatg[i1, j1]) - (phase - t * ωmatg[i2, j2])
            c20 = complex(a20 * Bm * cos(p20), a20 * Bm * sin(p20))
            v20 = complex(a20 * Am * cos(p20 - pi / 2), a20 * Am * sin(p20 - pi / 2))

            C22[iix2p, iiy2p] = C22[iix2p, iiy2p] + c22
            V22[iix2p, iiy2p] = V22[iix2p, iiy2p] + v22
            C20[iix2m, iiy2m] = C20[iix2m, iiy2m] + c20
            V20[iix2m, iiy2m] = V20[iix2m, iiy2m] + v20
        end
    end


    # Construct Fourier Components corresponding to second order

    F22 = 0.5 * prod(size(C22)) * ifftshift(ifftshift(C22, 1), 2)
    E22 = circshift(
        ifftshift(
            ifftshift(
                irfft(F22, (nx, ny)),
                1),
            2),
        (-1, -1))

    G22 = 0.5 * prod(size(V22)) * ifftshift(ifftshift(V22, 1), 2)
    P22 = circshift(
        ifftshift(
            ifftshift(
                irfft(G22, (nx, ny)),
                1),
            2),
        (-1, -1))

    F20 = 0.5 * prod(size(C20)) * ifftshift(ifftshift(C20, 1), 2)
    E20 = circshift(
        ifftshift(
            ifftshift(
                irfft(F20, (nx, ny)),
                1),
            2),
        (-1, -1))

    G20 = 0.5 * prod(size(V20)) * ifftshift(ifftshift(V20, 1), 2)
    P20 = reverse(circshift(
            ifftshift(
                ifftshift(
                    irfft(G20, (nx, ny)),
                    1),
                2),
            (-1, -1)),
        1)
    return E22, P22, E20, P20
end

end