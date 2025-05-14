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
    open_EP

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

    dirg = dirg + mcallister_mwd.(fmatg / fm, oi.twist_angle) # MWD as a function of freq
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

    kxmatg, kymatg, ampg_newwave_norm
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

function export_ow3d_init(η, ϕ, stime, oi::OW3DInput, dir)
    println("generating initial file...")
    file_name = "OceanWave3D_$(oi.nx)x$(oi.ny)_$(round(Int,100*oi.A))cm_rot$(round(Int,oi.twist_angle))_phase$(round(Int,oi.ϕ))_depth$(round(Int,oi.depth))_mwd$(round(Int,oi.mwd)).init"
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

function export_ow3d_inp(oi::OW3DInput, dir)
    println("generating input file...")
    file_name = "OceanWave3D_$(oi.nx)x$(oi.ny)_$(round(Int,100*oi.A))cm_rot$(round(Int,oi.twist_angle))_phase$(round(Int,oi.ϕ))_depth$(round(Int,oi.depth))_mwd$(round(Int,oi.mwd)).inp"
    Lx = (oi.nx - 1) * oi.dx
    Ly = (oi.ny - 1) * oi.dy
    fpath = joinpath(dir, file_name)
    open(fpath, "w") do file
        write(file, @sprintf "Tropical Cyclone focus wave H=%f nx=%d ny=%d dx=%f dy=%f depth=%f phase=%f Twist=%f MWD=%f" oi.A oi.nx oi.ny oi.dx oi.dy oi.depth oi.ϕ oi.twist_angle oi.mwd)

        write(file, @sprintf "A flat bottom, focused wave initial condition in 3D")
        write(file, @sprintf "-1 2                                  <- Initial condition (0=defined by funPressureTerm.f90, 1=NL standing wave, 2=shoaling on a smooth beach, 3=Whalin bar, ... see Initialization.f90:SetupInitialConditions); IncWaveType (0=none, 1=stream function, 2=linear regular or irregular waves)")
        write(file, @sprintf "%f %f %f %d %d 9 0 0 1 1 1 1   <- Lx Ly Lz Nx Ny Nz GridX GridY GridZ(0=even,1=clustering) GhostGrid (0=off,1=on)")
        write(file, @sprintf "3 3 3 1 1 1                          <- alpha, beta, gamma")
        write(file, @sprintf "%d %f 1 0 1                      <- Nsteps, dt, timeintegration scheme (1=RK4,2=lowstorage-RK45), CFL (if CFL/=0 then dt=CFL*dxmin/c, assume c=sqrt(g*hdeep)), RK4-ExtrapolationON/OFF")
        write(file, @sprintf "9.81                                 <- gravitational acceleration constant")
        write(file, @sprintf "1 3 0 55 1e-6 1e-6 1 V 1 1 20         <- solver (0:DC, 1:GMRES), GMRES Preconditioning (0=none (Matrix free,DIRECT),1=Linear LU(no elimination),2=Linear LU(ghostpoints eliminated),3=Multigrid (no elimination) ), Coarsening Strategy (0=Allan's, 1=Ole's), GMRESmaxiterations, relative tolerance (reltol), abstol, maxit, cyclet, pre-smoothings, post-smoothings, MGmaxgrids, DOF breakeven")
        write(file, @sprintf "0.05 1.00 1.84 2 0 0 1 6 32          <- Stream function solution parameters: H, h, L, T, WAVELorPER, uEorS, EorS, nsteps, maxiter ")
        write(file, @sprintf "5 20 1 3 <- StoreDataOnOff        <- StoreDataOnOff, formattype, (StoreDataOnOff=0 -> no output, StoreDataOnOff=+stride-> binary, StoreDataOnOff=-stride -> ascii every stride time steps.  formattype=0, binary; =1, unformatted) If formattype=20, then the line should read: StoreDataOnOff, iKinematics, formattype, nOutFiles; and nOutFiles lines should appear below defining  [xbeg, xend, xstride, ybeg, yend, ystride, tbeg, tend, tstride] for each file.")
        for kinematic in kinematics
            write(file, @sprintf "960 1080 1 240 250 1 1 1201 1  <- xbeg, xend, xstride, ybeg, yend, ystride, tbeg, tend, tstride")
        end
        write(file, @sprintf "1 0                                  <- 0/1=linear/nonlinear computations, Dynamic pressure term on/off")
        write(file, @sprintf "0 6 10 0.08 0.08 0.4                 <- SG-filtering on/off, filter half width, poly order")
        write(file, @sprintf "0 8. 3 Y 0.0                         <- relaxation zones on/off, transient time, no. zones. For each zone define on following lines: x1 x2 y1 y2 ftype(=relaxation function) param XorY WavegenONOFF Degrees(=IC rotation)")
        write(file, @sprintf "0 0                                  <- Damping pressure zone:  PDampingOnOff=0 (off), number of zones.  For each zone include the line: x1, x2, y1, y2 (bounds of the zone), gamma0 (dynamic FSBC), Gamma0 (kinematic FSBC), i_damperType (0=friction on the velocity, 1=friction on the potential).")
        write(file, @sprintf "0 2.0 2 0 0 1 0                      <- SWENSE on/off, ramp in time, wf direction (1:+x ; -1:-x ; 2:+y ; -2:-y ; >3: angle of the 3D wavefield), Reflexion of incident wf: West, East, North, South (0=off,1=on)")
        write(file, @sprintf "0                                    <- Curvilinear on/off ")
        write(file, @sprintf "33  8. 2. 80. 20. -1 -11 100. 50. run06.el 22.5 1.0 3.3 <- Irregular/regular waves:  i_spec, T_p, H_s, h0, kh_max, seed, seed2, x0, y0, (inc_wave_file or gamma_JONSWAP, if ispec=2 or 3), (inc_wave_file,beta,S,gamma_JONSWAP if ispec>=30).  For a random wave, the spectrum (i_spec=):  -1=>Monochromatic, 0=>P-M, 1=>JONSWAP, 2=>Read from a file, 3=>JONSWAP with input gamma; +- 3* means 3D at angle beta -30=>Monochromatic, 30=>P-M, 31=>JONSWAP, 32=>Not yet implemented 33=>JONSWAP with Normal spreading, 34=> JONSWAP with cos^S spreading.  ")
    end
end

function open_Kinematics(fpath; Nbits=32, compute_pressure=true)
    io_k = open(fpath)

    # Choose number of bits
    if Nbits == 32
        int_nbit = "int"
    elseif Nbits == 64
        int_nbit = "int64"
    else
        println("Illegal value $Nbits for Nbits; Nbits ∈ {32, 64}")
    end
    #
    # Read the data from the file
    # These read statements must correspond exactly to what appears in the
    # Fortran subroutine: <top dir>/src/IO/StoreKinematicData.f90
    #
    junk = fread(fid1, 1, int_nbit)
    xbeg = fread(fid1, 1, "int")
    xend = fread(fid1, 1, "int")
    xstride = fread(fid1, 1, "int")
    ybeg = fread(fid1, 1, "int")
    yend = fread(fid1, 1, "int")
    ystride = fread(fid1, 1, "int")
    tbeg = fread(fid1, 1, "int")
    tend = fread(fid1, 1, "int")
    tstride = fread(fid1, 1, "int")
    dt = fread(fid1, 1, "double") # Time step size
    nz = fread(fid1, 1, "int")
    junk = fread(fid1, 2, int_nbit) # Junk read statements are necessary for eol markers
    nx = floor((xend - xbeg) / xstride) + 1
    ny = floor((yend - ybeg) / ystride) + 1
    nt = floor((tend - tbeg) / tstride) + 1

    # A scratch vector for reading the data
    #
    tmp = zeros(nx * ny * max(nz, 5), 1)
    #
    # The x-y grid, the depth and bottom gradients for this slice of data
    #
    tmp(1:5*nx*ny) = fread(fid1, 5 * nx * ny, "double")
    junk = fread(fid1, 2, int_nbit)
    #
    x = zeros(nx, ny)
    x(:) = tmp(1:5:5*nx*ny)
    y = zeros(nx, ny)
    y(:) = tmp(2:5:5*nx*ny)
    h = zeros(nx, ny)
    h(:) = tmp(3:5:5*nx*ny)
    hx = zeros(nx, ny)
    hx(:) = tmp(4:5:5*nx*ny)
    hy = zeros(nx, ny)
    hy(:) = tmp(5:5:5*nx*ny)
    #
    # The sigma coordinate
    #
    for i = 1:nz
        sigma(i) = fread(fid1, 1, "double")
    end
    junk = fread(fid1, 2, int_nbit)
    #
    # Initialize arrays for the solution on this slice
    #
    eta = zeros(nt, nx, ny)
    etax = zeros(nt, nx, ny)
    etay = zeros(nt, nx, ny)
    phi = zeros(nt, nz, nx, ny)
    w = zeros(nt, nz, nx, ny)
    u = zeros(nt, nz, nx, ny)
    uz = zeros(nt, nz, nx, ny)
    v = zeros(nt, nz, nx, ny)
    vz = zeros(nt, nz, nx, ny)
    wz = zeros(nt, nz, nx, ny)
    t = [0:nt-1] * dt * tstride   # The time axis
    #
    # Read in the solution variables eta, gradeta, phi, u, v, w, dudz, dvdz.  
    #
    for it = 1:nt
        try
            tmp(1:nx*ny) = fread(fid1, nx * ny, "double")
            eta(it, :) = tmp(1:nx*ny)
            junk = fread(fid1, 2, int_nbit)
        catch
            warning(["Read failed at time step ", num2str(it)])
            break
        end
        #
        tmp(1:nx*ny) = fread(fid1, nx * ny, "double")
        etax(it, :) = tmp(1:nx*ny)
        junk = fread(fid1, 2, int_nbit)
        #
        tmp(1:nx*ny) = fread(fid1, nx * ny, "double")
        etay(it, :) = tmp(1:nx*ny)
        junk = fread(fid1, 2, int_nbit)
        #
        tmp = fread(fid1, nx * ny * nz, "double")
        phi(it, :) = tmp
        junk = fread(fid1, 2, int_nbit)
        #
        tmp = fread(fid1, nx * ny * nz, "double")
        u(it, :) = tmp
        junk = fread(fid1, 2, int_nbit)
        #
        tmp = fread(fid1, nx * ny * nz, "double")
        v(it, :) = tmp
        junk = fread(fid1, 2, int_nbit)
        #
        tmp = fread(fid1, nx * ny * nz, "double")
        w(it, :) = tmp
        junk = fread(fid1, 2, int_nbit)
        #
        [tmp, count] = fread(fid1, nx * ny * nz, "double")
        # Check for an incomplete run.
        if count == 0
            it = it - 1
            break
        end
        wz(it, :) = tmp
        junk = fread(fid1, 2, int_nbit)
        #
        [tmp, count] = fread(fid1, nx * ny * nz, "double")
        # Check for an incomplete run.
        if count == 0
            it = it - 1
            break
        end
        uz(it, :) = tmp
        junk = fread(fid1, 2, int_nbit)
        #
        [tmp, count] = fread(fid1, nx * ny * nz, "double")
        # Check for an incomplete run.
        if count == 0
            it = it - 1
            break
        end
        vz(it, :) = tmp
        junk = fread(fid1, 2, int_nbit)
    end
    display(["Read ',num2str(it),' data points out of ", num2str(nt)])

    ##
    if Pressure
        #
        # Compute the pressure and acceleration from the standard output
        # kinematics.  This is only done along the first slice in y for 3D
        # problems.
        #
        # Build the 4th-order even grid time differentiation matrix
        #
        alpha = 2
        r = 2 * alpha + 1
        c = BuildStencilEven(alpha, 1)
        Dt = spdiags([ones(nt, 1) * c(:, alpha + 1)'], [-alpha:alpha], nt, nt)
        for j = 1:alpha
            Dt(j, :) = 0
            Dt(j, 1:r) = c(:, j)'
            Dt(nt - j + 1, :) = 0
            Dt(nt - j + 1, nt-r+1:nt) = c(:, r - j + 1)'
        end
        Dt = Dt / dt
        #
        # Compute time-derivatives of eta, phi, and u and eta_tt
        #
        # ip=input("x point index to work with?");
        etat_m = zeros(nt, nx, ny)
        etatt_m = etat_m
        phit_m = zeros(nt, nz, nx, ny)
        p_m = phit_m
        ut_m = p_m
        for idy = 1:ny
            etat = zeros(nt, nx)
            etatt = etat
            phit = zeros(nt, nz, nx)
            p = phit
            ut = p
            ut1 = p
            for ip = 1:nx
                etat(:, ip) = Dt * eta(:, ip)
                etatt(:, ip) = Dt * etat(:, ip)
                #
                for j = 1:nz
                    phit(:, j, ip) = Dt * phi(:, j, ip) - w(:, j, ip) .* sigma(j) .* etat(:, ip)
                    p(:, j, ip) = -(phit(:, j, ip) + 1 / 2 * (u(:, j, ip) .^ 2 + w(:, j, ip) .^ 2))
                    ut(:, j, ip) = Dt * u(:, j, ip) - uz(:, j, ip) .* sigma(j) .* etat(:, ip)
                    ut1(:, j, ip) = Dt * u(:, j, ip)
                end
            end
            etat_m(:, :, idy) = etat
            etatt_m(:, :, idy) = etatt
            phit_m(:, :, :, idy) = phit
            p_m(:, :, :, idy) = p
            ut_m(:, :, :, idy) = ut
        end
    end
    it, eta, etat_m, etatt_m, etax, etay, phit_m, p_m, ut_m, u, v, w, uz, vz, wz, x, y, sigma, t
end
function BuildStencilEven(alpha, der)
    #
    # A function to compute finite-difference coefficients for the der^th 
    # derivative of a function on an evenly spaced grid.  2 alpha+1 sets of 
    # coefficients are returned where set 1,2,...,alpha are one-sided schemes 
    # at the left end, alpha+1 is the centered scheme, and alpha+2,...,rank 
    # are the one-sided schemes at the right end.  
    # 
    rank = 2 * alpha + 1
    # One-sided schemes for the left end-points.  
    for ip = 1:alpha
        for m = -ip+1:rank-ip
            for n = 1:rank
                mat(m + ip, n) = (m)^(n - 1) / factorial(n - 1)
            end
        end
        mat
        minv = inv(mat)
        fx(1:rank, ip) = minv(der + 1, 1:rank)'
    end
    # The centered scheme
    for m = -alpha:alpha
        for n = 1:rank
            mat(m + alpha + 1, n) = (m)^(n - 1) / factorial(n - 1)
        end
    end
    minv = inv(mat)
    fx(1:rank, alpha + 1) = minv(der + 1, 1:rank)'
    # Reflect the one-sided schemes from the left end.
    if mod(der, 2) == 0
        for ip = 1:alpha
            fx(1:rank, rank - ip + 1) = flipud(fx(1:rank, ip))
        end
    else
        for ip = 1:alpha
            fx(1:rank, rank - ip + 1) = -flipud(fx(1:rank, ip))
        end
    end
end
end
