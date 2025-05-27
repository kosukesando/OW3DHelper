module OW3DHelper
# using SpecialFunctions
using ProgressLogging
using TerminalLoggers
using LoopVectorization
using Printf
using Dates

export
    OW3DInput,
    calc_etaphi,
    JSpec,
    export_ow3d_init,
    open_EP,
    generate_init

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
    twist_type::String
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

    kxmatg, kymatg, ωmatg, ampg_newwave_norm
end

function calc_etaphi(oi::OW3DInput, t::Number)
    kxmatg, kymatg, ωmatg, ampg_newwave_norm = calc_k_omega_a(oi)
    η = calc_eta(oi, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm)
    ϕ = calc_phi(oi, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm, η)
    η, ϕ
end

function calc_etaphi(oi::OW3DInput, ts::Int, te::Int, x, y)
    kxmatg, kymatg, ωmatg, ampg_newwave_norm = calc_k_omega_a(oi)
    t = ts:te
    η = calc_eta(oi, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm, x, y)
    # ϕ = calc_phi(oi, kxmatg, kymatg, ωmatg, t, ampg_newwave_norm, x, y)
    ϕ = zeros(length(t))
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

function calc_eta(oi, kxmatg, kymatg, ωmatg, t_vec, ampg_newwave_norm, x::Float64, y::Float64)
    (nkx, nky) = size(ampg_newwave_norm)
    η = zeros(length(t_vec))
    for (i, t) in enumerate(t_vec)
        for kj = eachindex(1:nky), ki = eachindex(1:nkx)
            kx = kxmatg[ki, kj]
            ky = kymatg[ki, kj]
            ω = ωmatg[ki, kj]
            an = ampg_newwave_norm[ki, kj]
            phasei = kx * x + ky * y - ω * t + oi.ϕ
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
    ϕ = zeros(oi.nx, oi.ny)
    println("Calculating velocity potential at free surface")
    Threads.@threads for kj = eachindex(1:nky)
        for ki = eachindex(1:nkx)
            # for kj in eachindex(1:nky)
            kx = kxmatg[ki, kj]
            ky = kymatg[ki, kj]
            k = sqrt(kx^2 + ky^2)
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

function calc_etaphi_origin(oi::OW3DInput, ts::Int, te::Int)
    calc_etaphi(oi, ts, te, 0.0, 0.0)
end

function export_ow3d_init(η, ϕ, stime, oi::OW3DInput, dir; include_param=true)
    println("generating initial file...")
    if include_param
        file_name = "OceanWave3D_$(oi.nx)x$(oi.ny)_$(round(Int,100*oi.A))cm_rot$(round(Int,oi.twist_angle))_phase$(round(Int,oi.ϕ))_depth$(round(Int,oi.depth))_mwd$(round(Int,oi.mwd)).init"
    else
        file_name = "OceanWave3D.init"
    end
    Lx = (oi.nx - 1) * oi.dx
    Ly = (oi.ny - 1) * oi.dy
    fpath = joinpath(dir, file_name)
    open(fpath, "w") do file
        write(file, @sprintf "Tropical Cyclone focus wave H=%f nx=%d ny=%d dx=%f dy=%f depth=%f phase=%f Twist=%f MWD=%f" oi.A oi.nx oi.ny oi.dx oi.dy oi.depth oi.ϕ oi.twist_angle oi.mwd)
        write(file, @sprintf "\n%12e %12e %d %d %12e" Lx Ly oi.nx oi.ny stime)
        for ry = 1:oi.ny
            for rx = 1:oi.nx
                write(file, @sprintf "\n%12e %12e" η[rx, ry] ϕ[rx, ry])
            end
        end
    end
end

function generate_init(A::Float64, ϕ::Float64, k0::Float64, kmaxx::Float64, kmaxy::Float64, depth::Float64, dx::Float64, dy::Float64, nx::Int, ny::Int, spectrum::String, spreading_type::String, spreading_param::Float64, twist_angle::Float64, mwd::Float64, t::Float64; twist_type::String="mcallister", γ=3.3, dir=".")
    if spectrum == "JONSWAP"
        spec = JSpec(γ)
    end

    oi = OW3DInput(A, ϕ, k0, kmaxx, kmaxy, depth, dx, dy, nx, ny, spec, spreading_type, spreading_param, twist_angle, mwd, twist_type)
    eta, phi = calc_etaphi(oi, t)
    export_ow3d_init(eta, phi, 0, oi, dir; include_param=false)
end

function open_EP(fpath)
    io_ep = open(fpath)
    seek(io_ep, sizeof(Int32))
    Nx = read(io_ep, Int32)
    Ny = read(io_ep, Int32)
    seek(io_ep, sizeof(Int32) * 5)
    X = Array{Float64}(undef, Nx * Ny)
    Y = Array{Float64}(undef, Nx * Ny)
    read!(io_ep, X)
    read!(io_ep, Y)
    E = Array{Float64}(undef, Nx * Ny)
    P = Array{Float64}(undef, Nx * Ny)
    read!(io_ep, E)
    read!(io_ep, P)
    close(io_ep)
    Nx, Ny, X, Y, E, P
end

end
