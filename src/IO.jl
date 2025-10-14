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
    io = open(fpath)
    skip(io, sizeof(Int32))
    Nx = read(io, Int32)
    Ny = read(io, Int32)
    skip(io, sizeof(Int32) * 2)
    X = Array{Float64}(undef, Nx * Ny)
    Y = Array{Float64}(undef, Nx * Ny)
    read!(io, X)
    read!(io, Y)
    E = Array{Float64}(undef, Nx * Ny)
    P = Array{Float64}(undef, Nx * Ny)
    read!(io, E)
    read!(io, P)
    close(io)
    E = reshape(E, Nx, Ny)[2:end-1, 2:end-1]
    P = reshape(P, Nx, Ny)[2:end-1, 2:end-1]
    EPFile(Nx - 2, Ny - 2, reshape(X, Nx, Ny)[2:end-1, 2:end-1], reshape(Y, Nx, Ny)[2:end-1, 2:end-1], E, P)
end

function open_Kinematics(fpath)
    io = open(fpath)
    # This script reads the unformatted binary kinematics output file from the
    # OceanWave3D code.  
    #
    # Read the data from the file
    # These read statements must correspond exactly to what appears in the
    # Fortran subroutine: <top dir>/src/IO/StoreKinematicData.f90
    #
    skip(io, sizeof(Int32))
    xbeg = read(io, Int32) #
    xend = read(io, Int32) #
    xstride = read(io, Int32) #
    ybeg = read(io, Int32) #
    yend = read(io, Int32) #
    ystride = read(io, Int32) #
    tbeg = read(io, Int32) #
    tend = read(io, Int32) #
    tstride = read(io, Int32) #
    dt = read(io, Float64) # Time step size
    nz = read(io, Int32) #
    skip(io, 2 * sizeof(Int32))
    nx = floor(Int, (xend - xbeg) / xstride) + 1
    ny = floor(Int, (yend - ybeg) / ystride) + 1
    nt = floor(Int, (tend - tbeg) / tstride) + 1
    # The x-y grid, the depth and bottom gradients for this slice of data
    # tmp = Vector{Float64}(undef, max(nz, 5) * nx * ny)
    tmp = zeros(Float64, 5 * nx * ny)
    read!(io, tmp)
    skip(io, 2 * sizeof(Int32))
    x = reshape(tmp[1:5:5*nx*ny], nx, ny)
    y = reshape(tmp[2:5:5*nx*ny], nx, ny)
    h = reshape(tmp[3:5:5*nx*ny], nx, ny)
    hx = reshape(tmp[4:5:5*nx*ny], nx, ny)
    hy = reshape(tmp[5:5:5*nx*ny], nx, ny)
    sigma = zeros(nz)
    for i = 1:nz
        sigma[i] = read(io, Float64)
    end
    skip(io, 2 * sizeof(Int32))
    # Initialize arrays for the solution on this slice
    #
    eta = Array{Float64}(undef, nt, nx, ny)
    etax = Array{Float64}(undef, nt, nx, ny)
    etay = Array{Float64}(undef, nt, nx, ny)
    phi = Array{Float64}(undef, nt, nz, nx, ny)
    w = Array{Float64}(undef, nt, nz, nx, ny)
    u = Array{Float64}(undef, nt, nz, nx, ny)
    uz = Array{Float64}(undef, nt, nz, nx, ny)
    v = Array{Float64}(undef, nt, nz, nx, ny)
    vz = Array{Float64}(undef, nt, nz, nx, ny)
    wz = Array{Float64}(undef, nt, nz, nx, ny)
    t = collect(0:nt-1) * dt * tstride   # The time axis
    #
    # Read in the solution variables eta, gradeta, phi, u, v, w, dudz, dvdz.  
    #
    for it = 1:nt
        read!(io, @views eta[it, :, :])
        skip(io, 2 * sizeof(Int32))
        read!(io, @views etax[it, :, :])
        skip(io, 2 * sizeof(Int32))
        read!(io, @views etay[it, :, :])
        skip(io, 2 * sizeof(Int32))
        read!(io, @views phi[it, :, :, :])
        skip(io, 2 * sizeof(Int32))
        read!(io, @views u[it, :, :, :])
        skip(io, 2 * sizeof(Int32))
        read!(io, @views v[it, :, :, :])
        skip(io, 2 * sizeof(Int32))
        read!(io, @views w[it, :, :, :])
        skip(io, 2 * sizeof(Int32))
        read!(io, @views wz[it, :, :, :])
        skip(io, 2 * sizeof(Int32))
        read!(io, @views uz[it, :, :, :])
        skip(io, 2 * sizeof(Int32))
        read!(io, @views vz[it, :, :, :])
        skip(io, 2 * sizeof(Int32))
    end
    @info("Read ", nt, "data points")
    return KinematicsFile(
        xbeg,
        xend,
        xstride,
        ybeg,
        yend,
        ystride,
        tbeg,
        tend,
        tstride,
        dt,
        nz,
        nx,
        ny,
        nt,
        sigma,
        t,
        x,
        y,
        h,
        hx,
        hy,
        eta,
        etax,
        etay,
        phi,
        u,
        uz,
        v,
        vz,
        w,
        wz,
    )
end

function calc_pressure(kf::KinematicsFile)
    nx = kf.nx
    ny = kf.ny
    nz = kf.nz
    nt = kf.ny
    dt = kf.dt
    sigma = kf.sigma
    eta = kf.eta
    phi = kf.phi
    u = kf.u
    uz = kf.uz
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
        mat = zeros()
        fx = zeros()
        for ip = 1:alpha
            for m = -ip+1:rank-ip
                for n = 1:rank
                    mat[m+ip, n] = (m)^(n - 1) / factorial(n - 1)
                end
            end
            minv = inv(mat)
            fx[1:rank, ip] = minv(der + 1, 1:rank)'
        end
        # The centered scheme
        for m = -alpha:alpha
            for n = 1:rank
                mat[m+alpha+1, n] = (m)^(n - 1) / factorial(n - 1)
            end
        end
        minv = inv(mat)
        fx[1:rank, alpha+1] = minv(der + 1, 1:rank)'
        # Reflect the one-sided schemes from the left end.
        if mod(der, 2) == 0
            for ip = 1:alpha
                fx[1:rank, rank-ip+1] = reverse(fx[1:rank, ip], dims=1)
            end
        else
            for ip = 1:alpha
                fx[1:rank, rank-ip+1] = -reverse(fx[1:rank, ip], dims=1)
            end
        end
    end
    if calc_pressure
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
        Dt = spdiagm([ones(nt, 1) * c[:, alpha+1]'], [-alpha:alpha], nt, nt)
        for j = 1:alpha
            Dt[j, :] = 0
            Dt[j, 1:r] = c[:, j]'
            Dt[nt-j+1, :] = 0
            Dt[nt-j+1, nt-r+1:nt] = c[:, r-j+1]'
        end
        Dt = Dt / dt
        #
        # Compute time-derivatives of eta, phi, and u and eta_tt
        #
        # ip=input("x point index to work with?");
        etat_m = Array{Float64}(undef, nt, nx, ny)
        etatt_m = etat_m
        phit_m = Array{Float64}(undef, nt, nz, nx, ny)
        p_m = phit_m
        ut_m = p_m
        for idy = 1:ny
            etat = Array{Float64}(undef, nt, nx)
            etatt = etat
            phit = Array{Float64}(undef, nt, nz, nx)
            p = phit
            ut = p
            ut1 = p
            for ip = 1:nx
                etat[:, ip] = Dt * eta[:, ip]
                etatt[:, ip] = Dt * etat[:, ip]
                #
                for j = 1:nz
                    phit[:, j, ip] = Dt * phi[:, j, ip] - w[:, j, ip] .* sigma[j] .* etat[:, ip]
                    p[:, j, ip] = -(phit[:, j, ip] + 1 / 2 * (u[:, j, ip] .^ 2 + w[:, j, ip] .^ 2))
                    ut[:, j, ip] = Dt * u[:, j, ip] - uz[:, j, ip] .* sigma[j] .* etat[:, ip]
                    ut1[:, j, ip] = Dt * u[:, j, ip]
                end
            end
            etat_m[:, :, idy] = etat
            etatt_m[:, :, idy] = etatt
            phit_m[:, :, :, idy] = phit
            p_m[:, :, :, idy] = p
            ut_m[:, :, :, idy] = ut
        end
    end
end

function getEPlist(casename::String, twist::Int, phase::Int; basedir=".")
    epfiledirs = []
    for i in 1:5
        dir = joinpath(basedir, @sprintf("%s/%03ddeg/%03d/%d", casename, twist, phase, i))
        if isdir(dir)
            if i == 1
                append!(epfiledirs, filter(!contains("EP_99999.bin"), filter(contains(r"EP_*"), readdir(realpath(dir), join=true))))
            else
                ep_head = joinpath(dir, "EP_00000.bin")
                if !isfile(ep_head)
                    break
                end
                ep1 = open_EP(epfiledirs[end])
                ep2 = open_EP(ep_head)
                if ep1.η == ep2.η
                    append!(epfiledirs, filter(!contains("EP_00000.bin"), filter(!contains("EP_99999.bin"), filter(contains(r"EP_*"), readdir(realpath(dir), join=true)))))
                else
                    append!(epfiledirs, filter(!contains("EP_99999.bin"), filter(contains(r"EP_*"), readdir(realpath(dir), join=true))))
                end
            end
        end
    end
    epfiledirs
end

getEPlist(s; basedir=".") = getEPlist(s.casename, s.twist, s.phase)
