module OW3DHelper
# using SpecialFunctions
using ProgressLogging
using TerminalLoggers
using LoopVectorization

export
    OW3DInput,
    calc_etaphi,
    JSpec,
    calc_eta_origin,

const g::Float64 = 9.81

abstract type AbstractSpecType end

struct JSpec <: AbstractSpecType
    γ::Float64
end

struct GSpec <: AbstractSpecType
    num_anglc
    num_specc
end

struct OW3DInput{T}
    A::Float64
    ϕ::Float64
    k0::Float64
    kmaxx::Float64
    kmaxy::Float64
    depth::Float64
    dx::Float64
    dy::Float64
    nx::Int
    ny::Int
    # stime::Float64
    spec::T
    spreading_type::String
    spreading_param::Float64
    twist_angle::Float64
    mwd::Float64
end

function mcallister_mwd(foverfp::Float64, skewd::Float64)
    # Based on parameterisation in McAllister(2019)
    if foverfp <= 1
        theta = 0.0
    else
        theta = (tanh(2 * foverfp - 3.5) / tanh(1.5) + 1) * skewd / 2
    end
    return theta
end

function calc_etaphi(oi::OW3DInput, t::Number)
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
        if oi.spreading_type == "T"
            Sdir = exp(-0.5 * (dirg .^ 2) ./ (oi.spreading_param^2))
        end
        if oi.spreading_type == "E"

            s = zeros(1, len)
            tm1 = zeros(1, len)
            tm2 = zeros(1, len)
            BM = zeros(1, len)

            IL = fmatg .< fm
            IG = fmatg .>= fm

            s[IL] .= 0.0 .+ 1.0 * ((fmatg[IL] ./ fm) .^ (-7.929))
            s[IG] .= 0.0 .- 1.0 * ((fmatg[IG] ./ fm) .^ (-2))

            tm1[IL] .+= 14.93 / 2
            tm2[IL] .-= 14.93 / 2

            tm1[IG] .+= (1 / 2) * (exp.(5.453 .- 2.750 .* ((fmatg[IG] ./ fm) .^ (-1))))
            tm2[IG] .-= (1 / 2) * (exp.(5.453 .- 2.750 .* ((fmatg[IG] ./ fm) .^ (-1))))

            E = (1 / (sqrt(8 * pi))) * (1 ./ s)

            for kk = -10:1:10
                BM1 = exp.((-0.5) * ((((pi / 180) * dirg .- (pi / 180) * tm1 .- 2 * pi * kk) ./ ((pi / 180) * s)) .^ 2))
                BM2 = exp.((-0.5) * ((((pi / 180) * dirg .- (pi / 180) * tm2 .- 2 * pi * kk) ./ ((pi / 180) * s)) .^ 2))
                BM .+= BM1 + BM2
            end

            Sdir = E .* BM # Ewans bi-modal spreading function

        end

        ampg_newwave = 1 ./ kmatg .* Spec .* Sdir #amplitude of component - S(kx,ky)
        replace!(ampg_newwave, NaN => 0)
    else
        println("spectrum not recognised")
    end

    dirg = dirg + mcallister_mwd.(fmatg / fm, oi.twist_angle) # MWD as a function of freq
    dirg = dirg .+ oi.mwd
    kxmatg = kmatg .* cos.(dirg * (pi / 180)) #Component of k in x direction
    kymatg = kmatg .* sin.(dirg * (pi / 180)) #Component of k in y direction
    # wtmatg = wmatg * oi.stime       #omega * t
    wtmatg = wmatg * t       #omega * t

    ksmx = oi.kmaxx - 1.00 * oi.k0 # start smoothing kx truncation
    ksmxl = 0.5 * oi.k0 # start smoothing kx lower-bound
    Ismx = abs.(kxmatg) .>= ksmx .&& abs.(kxmatg) .<= oi.kmaxx
    Ismxl = abs.(kxmatg) .>= kminx .&& abs.(kxmatg) .<= ksmxl
    # ampg_newwave[Ismx] .= ampg_newwave[Ismx] .* (0.5 * (cos.(((abs.(kxmatg[Ismx]) .- ksmx) ./ (oi.kmaxx .- ksmx)) * (pi)) .+ 1)) .^ (1)
    # ampg_newwave[Ismxl] .= ampg_newwave[Ismxl] .* (0.5 * (cos.(((abs.(kxmatg[Ismxl]) .- ksmxl) ./ (kminx .- ksmxl)) * (pi)) .+ 1)) .^ (1)
    # ampg_newwave[abs.(kxmatg).>=oi.kmaxx] .= 0
    # ampg_newwave[abs.(kxmatg).<=kminx] .= 0

    ksmy = oi.kmaxy - 1.00 * oi.k0 # start smoothing ky truncation
    ksmyl = 0.5 * oi.k0 # start smoothing ky lower-bound
    Ismy = abs.(kymatg) .>= ksmy .&& abs.(kymatg) .<= oi.kmaxy
    Ismyl = abs.(kymatg) .>= kminy .&& abs.(kymatg) .<= ksmyl
    # ampg_newwave[Ismy] .= ampg_newwave[Ismy] .* (0.5 * (cos.(((abs.(kymatg[Ismy]) .- ksmy) ./ (oi.kmaxy .- ksmy)) * (pi)) .+ 1)) .^ (1)
    # ampg_newwave[Ismyl] .= ampg_newwave[Ismyl] .* (0.5 * (cos.(((abs.(kymatg[Ismyl]) .- ksmyl) ./ (kminy .- ksmyl)) * (pi)) .+ 1)) .^ (1)
    # ampg_newwave[abs.(kymatg).>=oi.kmaxy] .= 0
    # ampg_newwave[abs.(kymatg).<=kminy] .= 0

    ampg_newwave_norm = oi.A * (ampg_newwave ./ sum(ampg_newwave)) #  normalise to control size of event
    # return ampg_newwave_norm
    (nkx, nky) = size(ampg_newwave_norm)
    # Generate full-domain
    xvec = oi.dx * (-(oi.nx - 1)/2:1:(oi.nx-1)/2)
    yvec = oi.dy * (-(oi.ny - 1)/2:1:(oi.ny-1)/2)

    # η = zeros(oi.nx, oi.ny)
    # ϕ = zeros(oi.nx, oi.ny)
    η = calc_eta(oi, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm)
    # calc_eta!(η, oi, kxmatg, kymatg, wtmatg, ampg_newwave_norm)
    ϕ = calc_phi(oi, kmatg, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm, η)
    # ϕ = zeros(oi.nx, oi.ny)
    η, ϕ
end

function calc_eta(oi, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm)
    (nkx, nky) = size(ampg_newwave_norm)
    η = zeros(oi.nx, oi.ny)
    # Generate full-domain
    xvec = oi.dx * (-(oi.nx - 1)/2:1:(oi.nx-1)/2)
    yvec = oi.dy * (-(oi.ny - 1)/2:1:(oi.ny-1)/2)
    # Calculate linear free surface
    println("Calculating linear free surface")
    Threads.@threads for kj = eachindex(1:nky)
        for ki = eachindex(1:nkx)
            kx = kxmatg[ki, kj]
            ky = kymatg[ki, kj]
            ω = ωmatg[ki, kj]
            an = ampg_newwave_norm[ki, kj]
            for yi = eachindex(yvec), xi = eachindex(xvec)
                x = xvec[xi]
                y = yvec[yi]
                phasei = kx * x + ky * y - ω * t + oi.ϕ
                etacomp = an * cos(phasei)
                # if !isnan(etacomp)
                η[xi, yi] += etacomp
                # end
            end
        end
    end
    return η
end

function _calc_eta_single(oi, kxmatg, kymatg, wtmatg_arr, ampg_newwave_norm, x::Float64, y::Float64)
    (nkx, nky) = size(ampg_newwave_norm)
    # Generate full-domain
    # Calculate linear free surface
    println("Calculating linear free surface")
    η = zeros(length(wtmatg_arr))
    for (i, wtmatg) in enumerate(wtmatg_arr)
        for kj = eachindex(1:nky), ki = eachindex(1:nkx)
            kx = kxmatg[ki, kj]
            ky = kymatg[ki, kj]
            wt = wtmatg[ki, kj]
            an = ampg_newwave_norm[ki, kj]
            phasei = kx * x + ky * y - wt + oi.ϕ
            etacomp = an * cos(phasei)
            η[i] += etacomp
        end
    end
    η
end

function calc_phi(oi, kmatg, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm, η)
    (nkx, nky) = size(ampg_newwave_norm)
    # Generate full-domain
    xvec = oi.dx * (-(oi.nx - 1)/2:1:(oi.nx-1)/2)
    yvec = oi.dy * (-(oi.ny - 1)/2:1:(oi.ny-1)/2)
    ϕ = zeros(oi.nx, oi.ny)
    println("Calculating velocity potential at free surface")
    Threads.@threads for kj = eachindex(1:nky)
        for ki = eachindex(1:nkx)
            # for kj in eachindex(1:nky)
            k = kmatg[ki, kj]
            kx = kxmatg[ki, kj]
            ky = kymatg[ki, kj]
            ω = ωmatg[ki, kj]
            an = ampg_newwave_norm[ki, kj]
            for xi in eachindex(xvec)
                for yi in eachindex(yvec)
                    x = xvec[xi]
                    y = yvec[yi]
                    phasei = kx * x + ky * y - ω * t + oi.ϕ
                    phicomp = (an / (ω + 0.000000001)) * g * ((cosh(k * (η[xi, yi] + oi.depth)) / cosh(k * oi.depth)) * sin(phasei))
                    @views(ϕ[xi, yi] += phicomp)
                end
            end
        end
    end
    ϕ
end

function calc_eta_origin(oi::OW3DInput, ts::Int, te::Int)
    calc_eta_single(oi, ts, te, 0.0, 0.0)
end

function calc_eta_single(oi::OW3DInput, ts::Int, te::Int, x, y)
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
        # Distribution valid for S(w) & S(f) since scaled with normalisation

        if oi.spreading_type == "W"
            Sdir = exp.(-0.5 * (dirg .^ 2) ./ (oi.spreading_param^2)) # A wrapped normal spreading function with parameter sigma
        end
        if oi.spreading_type == "T"
            Sdir = exp(-0.5 * (dirg .^ 2) ./ (oi.spreading_param^2))
        end
        if oi.spreading_type == "E"

            s = zeros(1, len)
            tm1 = zeros(1, len)
            tm2 = zeros(1, len)
            BM = zeros(1, len)

            IL = fmatg .< fm
            IG = fmatg .>= fm

            s[IL] .= 0.0 .+ 1.0 * ((fmatg[IL] ./ fm) .^ (-7.929))
            s[IG] .= 0.0 .- 1.0 * ((fmatg[IG] ./ fm) .^ (-2))

            tm1[IL] .+= 14.93 / 2
            tm2[IL] .-= 14.93 / 2

            tm1[IG] .+= (1 / 2) * (exp.(5.453 .- 2.750 .* ((fmatg[IG] ./ fm) .^ (-1))))
            tm2[IG] .-= (1 / 2) * (exp.(5.453 .- 2.750 .* ((fmatg[IG] ./ fm) .^ (-1))))

            E = (1 / (sqrt(8 * pi))) * (1 ./ s)

            for kk = -10:1:10
                BM1 = exp.((-0.5) * ((((pi / 180) * dirg .- (pi / 180) * tm1 .- 2 * pi * kk) ./ ((pi / 180) * s)) .^ 2))
                BM2 = exp.((-0.5) * ((((pi / 180) * dirg .- (pi / 180) * tm2 .- 2 * pi * kk) ./ ((pi / 180) * s)) .^ 2))
                BM .+= BM1 + BM2
            end

            Sdir = E .* BM # Ewans bi-modal spreading function

        end

        ampg_newwave = 1 ./ kmatg .* Spec .* Sdir #amplitude of component - S(kx,ky)
        replace!(ampg_newwave, NaN => 0)
    else
        println("spectrum not recognised")
    end

    dirg = dirg + mcallister_mwd.(fmatg / fm, oi.twist_angle) # MWD as a function of freq
    dirg = dirg .+ oi.mwd
    kxmatg = kmatg .* cos.(dirg * (pi / 180)) #Component of k in x direction
    kymatg = kmatg .* sin.(dirg * (pi / 180)) #Component of k in y direction
    # wtmatg = wmatg * oi.stime       #omega * t

    # ksmx = oi.kmaxx - 1.00 * oi.k0 # start smoothing kx truncation
    # ksmxl = 0.5 * oi.k0 # start smoothing kx lower-bound
    # Ismx = abs.(kxmatg) .>= ksmx .&& abs.(kxmatg) .<= oi.kmaxx
    # Ismxl = abs.(kxmatg) .>= kminx .&& abs.(kxmatg) .<= ksmxl
    # ampg_newwave[Ismx] .= ampg_newwave[Ismx] .* (0.5 * (cos.(((abs.(kxmatg[Ismx]) .- ksmx) ./ (oi.kmaxx .- ksmx)) * (pi)) .+ 1)) .^ (1)
    # ampg_newwave[Ismxl] .= ampg_newwave[Ismxl] .* (0.5 * (cos.(((abs.(kxmatg[Ismxl]) .- ksmxl) ./ (kminx .- ksmxl)) * (pi)) .+ 1)) .^ (1)
    # ampg_newwave[abs.(kxmatg).>=oi.kmaxx] .= 0
    # ampg_newwave[abs.(kxmatg).<=kminx] .= 0

    # ksmy = oi.kmaxy - 1.00 * oi.k0 # start smoothing ky truncation
    # ksmyl = 0.5 * oi.k0 # start smoothing ky lower-bound
    # Ismy = abs.(kymatg) .>= ksmy .&& abs.(kymatg) .<= oi.kmaxy
    # Ismyl = abs.(kymatg) .>= kminy .&& abs.(kymatg) .<= ksmyl
    # ampg_newwave[Ismy] .= ampg_newwave[Ismy] .* (0.5 * (cos.(((abs.(kymatg[Ismy]) .- ksmy) ./ (oi.kmaxy .- ksmy)) * (pi)) .+ 1)) .^ (1)
    # ampg_newwave[Ismyl] .= ampg_newwave[Ismyl] .* (0.5 * (cos.(((abs.(kymatg[Ismyl]) .- ksmyl) ./ (kminy .- ksmyl)) * (pi)) .+ 1)) .^ (1)
    # ampg_newwave[abs.(kymatg).>=oi.kmaxy] .= 0
    # ampg_newwave[abs.(kymatg).<=kminy] .= 0

    ampg_newwave_norm = oi.A * (ampg_newwave ./ sum(ampg_newwave)) #  normalise to control size of event
    # return ampg_newwave_norm
    (nkx, nky) = size(ampg_newwave_norm)
    # Generate full-domain

    t = ts:te
    wtmatg_arr = [ωmatg * _t for _t in t]
    η = _calc_eta_single(oi, kxmatg, kymatg, wtmatg_arr, ampg_newwave_norm, x, y)
    # ϕ = calc_phi(oi, kmatg, kxmatg, kymatg, wtmatg, ampg_newwave_norm, η)
    η
end
end
