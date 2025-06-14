module OW3DHelper
# using SpecialFunctions
using Printf
using OhMyThreads: @tasks, @local, @set

export
    OW3DInput,
    EPFile,
    calc_etaphi,
    JSpec,
    KinematicSetting,
    export_ow3d_init,
    export_ow3d_inp,
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

struct KinematicSetting
    xbeg::Int
    xend::Int
    xstride::Int
    ybeg::Int
    yend::Int
    ystride::Int
    tbeg::Int
    tend::Int
    tstride::Int
end

struct EPFile
    nx::Int
    ny::Int
    x::Array{Float64}
    y::Array{Float64}
    η::Matrix{Float64}
    ϕ::Matrix{Float64}
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
        write(file, @sprintf "Tropical Cyclone focus wave H=%0.3f nx=%d ny=%d dx=%0.3f dy=%0.3f depth=%0.3f phase=%0.3f Twist=%0.3f MWD=%0.3f" oi.A oi.nx oi.ny oi.dx oi.dy oi.depth oi.ϕ oi.twist_angle oi.mwd)
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

function export_ow3d_inp(oi::OW3DInput, nt, kinematics, dir; include_param=true, nz=9, di_kine=5, dt=0.2)
    println("generating input file...")
    if include_param
        file_name = "OceanWave3D_$(oi.nx)x$(oi.ny)_$(round(Int,100*oi.A))cm_rot$(round(Int,oi.twist_angle))_phase$(round(Int,oi.ϕ))_depth$(round(Int,oi.depth))_mwd$(round(Int,oi.mwd)).inp"
    else
        file_name = "OceanWave3D.inp"
    end
    Lx = (oi.nx - 1) * oi.dx
    Ly = (oi.ny - 1) * oi.dy
    fpath = joinpath(dir, file_name)
    open(fpath, "w") do file
        write(file, "A flat bottom, focused wave initial condition in 3D\n")
        write(
            file,
            "-1 2 0.5                                 <- Initial condition (0=defined by funPressureTerm.f90, 1=NL standing wave, 2=shoaling on a smooth beach, 3=Whalin bar, ... see Initialization.f90:SetupInitialConditions); IncWaveType (0=none, 1=stream function, 2=linear regular or irregular waves)\n"
        )
        write(file, @sprintf "%0.3f %0.3f %0.3f %d %d %d 0 0 1 1 1 1    <- Lx Ly Lz Nx Ny Nz GridX GridY GridZ(0=even,1=clustering) GhostGrid (0=off,1=on)\n" Lx Ly oi.depth oi.nx oi.ny nz)
        write(
            file,
            "3 3 3 1 1 1                           <- alpha, beta, gamma\n"
        )
        write(file, @sprintf "%d %0.3f 1 0 1                       <- Nsteps, dt, timeintegration scheme (1=RK4,2=lowstorage-RK45), CFL (if CFL/=0 then dt=CFL*dxmin/c, assume c=sqrt(g*hdeep)), RK4-ExtrapolationON/OFF\n" nt dt)
        write(
            file,
            "9.81                                  <- gravitational acceleration constant\n"
        )
        write(
            file,
            "1 3 0 55 1e-6 1e-6 1 V 1 1 20          <- solver (0:DC, 1:GMRES), GMRES Preconditioning (0=none (Matrix free,DIRECT),1=Linear LU(no elimination),2=Linear LU(ghostpoints eliminated),3=Multigrid (no elimination) ), Coarsening Strategy (0=Allan's, 1=Ole's), GMRESmaxiterations, relative tolerance (reltol), abstol, maxit, cyclet, pre-smoothings, post-smoothings, MGmaxgrids, DOF breakeven\n"
        )
        write(
            file,
            "0.05 1.00 1.84 2 0 0 1 6 32           <- Stream function solution parameters: H, h, L, T, WAVELorPER, uEorS, EorS, nsteps, maxiter \n"
        )
        write(file, @sprintf "%d 20 1 %d         <- StoreDataOnOff, formattype, (StoreDataOnOff=0 -> no output, StoreDataOnOff=+stride-> binary, StoreDataOnOff=-stride -> ascii every stride time steps.  formattype=0, binary; =1, unformatted) If formattype=20, then the line should read: StoreDataOnOff, iKinematics, formattype, nOutFiles; and nOutFiles lines should appear below defining  [xbeg, xend, xstride, ybeg, yend, ystride, tbeg, tend, tstride] for each file.\n" di_kine length(kinematics))
        for kinematic in kinematics
            write(file, @sprintf "%d %d %d %d %d %d %d %d %d   <- xbeg, xend, xstride, ybeg, yend, ystride, tbeg, tend, tstride\n" kinematic.xbeg kinematic.xend kinematic.xstride kinematic.ybeg kinematic.yend kinematic.ystride kinematic.tbeg kinematic.tend kinematic.tstride)
        end
        write(file, " \n")
        write(
            file,
            "1 0                                   <- 0/1=linear/nonlinear computations, Dynamic pressure term on/off\n"
        )
        write(
            file,
            "0 6 10 0.08 0.08 0.4                  <- SG-filtering on/off, filter half width, poly order\n"
        )
        write(
            file,
            "0 8. 3 Y 0.0                          <- relaxation zones on/off, transient time, no. zones. For each zone define on following lines: x1 x2 y1 y2 ftype(=relaxation function) param XorY WavegenONOFF Degrees(=IC rotation)\n"
        )
        write(
            file,
            "0 0                                   <- Damping pressure zone:  PDampingOnOff=0 (off), number of zones.  For each zone include the line: x1, x2, y1, y2 (bounds of the zone), gamma0 (dynamic FSBC), Gamma0 (kinematic FSBC), i_damperType (0=friction on the velocity, 1=friction on the potential).\n"
        )
        write(
            file,
            "0 2.0 2 0 0 1 0                       <- SWENSE on/off, ramp in time, wf direction (1:+x ; -1:-x ; 2:+y ; -2:-y ; >3: angle of the 3D wavefield), Reflexion of incident wf: West, East, North, South (0=off,1=on)\n"
        )
        write(
            file,
            "0                                     <- Curvilinear on/off \n"
        )
        write(
            file,
            "33  8. 2. 80. 20. -1 -11 100. 50. run06.el 22.5 1.0 3.3  <- Irregular/regular waves:  i_spec, T_p, H_s, h0, kh_max, seed, seed2, x0, y0, (inc_wave_file or gamma_JONSWAP, if ispec=2 or 3), (inc_wave_file,beta,S,gamma_JONSWAP if ispec>=30).  For a random wave, the spectrum (i_spec=):  -1=>Monochromatic, 0=>P-M, 1=>JONSWAP, 2=>Read from a file, 3=>JONSWAP with input gamma; +- 3* means 3D at angle beta -30=>Monochromatic, 30=>P-M, 31=>JONSWAP, 32=>Not yet implemented 33=>JONSWAP with Normal spreading, 34=> JONSWAP with cos^S spreading.  \n"
        )
    end
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
    E = reshape(E, Nx, Ny)[2:end-1, 2:end-1]
    P = reshape(P, Nx, Ny)[2:end-1, 2:end-1]
    EPFile(Nx - 2, Ny - 2, reshape(X, Nx, Ny)[2:end-1, 2:end-1], reshape(Y, Nx, Ny)[2:end-1, 2:end-1], E, P)
end

end
