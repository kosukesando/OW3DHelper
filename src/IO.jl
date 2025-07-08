module IO

using Printf

export
    export_ow3d_init,
    export_ow3d_inp,
    open_EP,
    generate_init

include("./OW3DHelper.jl")
using .OW3DHelper

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