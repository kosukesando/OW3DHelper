function export_nc_ep(s; basedir=".", force=false)
    fname = joinpath(basedir, @sprintf("%s/%03ddeg/ep.nc", s.casename, s.twist))
    phases = [0, 90, 180, 270]
    epl = [getEPlist(s.casename, s.twist, p) for p in phases]
    if any(length.(epl) .< s.N) && !force
        @warn "Not enough EP files!"
        return
    end
    if isfile(fname)
        @warn "NC file exists!"
        return
    end

    ds = NCDataset(fname, "c")
    defDim(ds, "x", s.nx)
    defDim(ds, "y", s.ny)
    defDim(ds, "t", s.N)
    ds.attrib["casename"] = s.casename
    ds.attrib["twist"] = s.twist
    ds.attrib["twist_model"] = s.twist_model
    defVar(ds, "Lx", s.Lx, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "Ly", s.Ly, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "dx", s.dx, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "dy", s.dy, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "dt", s.dt, (), attrib=OrderedDict("units" => "s",))
    defVar(ds, "amplitude", s.amp, (), attrib=OrderedDict("units" => "m",))

    for i in 1:4
        η = zeros(s.N, s.nx, s.ny)
        ϕ = zeros(s.N, s.nx, s.ny)
        Threads.@threads for t in 1:s.N
            local ep = open_EP(epl[i][t])
            @views η[t, :, :] = ep.η
            @views ϕ[t, :, :] = ep.ϕ
        end

        defVar(ds, @sprintf("eta%03d", phases[i]), η, ("t", "x", "y"), attrib=OrderedDict(
            "units" => "m",
        ))
        defVar(ds, @sprintf("phi%03d", phases[i]), ϕ, ("t", "x", "y"), attrib=OrderedDict(
            "units" => "m²/s",
        ))
    end

    close(ds)
    @info @sprintf("Write complete for ep.nc %s(%ddeg twist)", s.casename, s.twist)
end

function export_nc_hilbert(s; basedir=".")
    phases = [0, 90, 180, 270]
    dir = joinpath(basedir, @sprintf("%s/%03ddeg", s.casename, s.twist))
    fname = joinpath(dir, "hilbert.nc")
    ds = NCDataset(fname, "c")
    defDim(ds, "x", s.nx)
    defDim(ds, "y", s.ny)
    defDim(ds, "t", s.N)
    ds.attrib["title"] = "hilbert transform"
    ds.attrib["casename"] = s.casename
    ds.attrib["twist"] = s.twist
    ds.attrib["twist_model"] = s.twist_model
    defVar(ds, "Lx", s.Lx, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "Ly", s.Ly, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "dx", s.dx, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "dy", s.dy, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "dt", s.dt, (), attrib=OrderedDict("units" => "s",))
    defVar(ds, "amplitude", s.amp, (), attrib=OrderedDict("units" => "m",))
    f_input = joinpath(dir, "ep.nc")
    if !isfile(f_input)
        @warn "ep.nc does not exist"
        return
    end
    ds_ep = NCDataset(f_input)
    for i in 1:4
        η = Array(ds_ep[@sprintf("eta%03d", phases[i])])
        ϕ = Array(ds_ep[@sprintf("phi%03d", phases[i])])
        η_hilbert = imag.(mapslices(hilbert, η; dims=1))
        ϕ_hilbert = imag.(mapslices(hilbert, ϕ; dims=1))

        defVar(ds, @sprintf("eta%03d", phases[i]), η_hilbert, ("t", "x", "y"), attrib=OrderedDict(
            "units" => "m",
        ))
        defVar(ds, @sprintf("phi%03d", phases[i]), ϕ_hilbert, ("t", "x", "y"), attrib=OrderedDict(
            "units" => "m²/s",
        ))
    end
    close(ds)
    @info @sprintf("Write complete for hilbert.nc: %s(%ddeg twist, phase=%d)", s.casename, s.twist)
end

function calc_argmax(s)
    epl = getEPlist(s)
    emax = 0
    emaxloc = nothing
    emaxt = 0
    for (i, ep) in enumerate(epl[1:s.N])
        e = open_EP(ep).η
        if maximum(e) > emax
            emaxloc = argmax(e)
            emax = maximum(e)
            emaxt = i
        end
    end
    (emaxt, Tuple(emaxloc)...)
end


function export_nc_4phase(s; basedir=".")
    dir = joinpath(basedir, @sprintf("%s/%03ddeg", s.casename, s.twist))
    f_input = joinpath(dir, "ep.nc")
    fh_input = joinpath(dir, "hilbert.nc")
    if !isfile(f_input)
        @warn "ep.nc does not exist"
        return
    end
    if !isfile(fh_input)
        @warn "hilbert.nc does not exist"
        return
    end
    ds_ep = NCDataset(f_input)
    ds_eph = NCDataset(fh_input)

    fname = joinpath(dir, "4p.nc")
    ds = NCDataset(fname, "c")
    defVar(ds, "dt", s.dt, (), attrib=OrderedDict("units" => "s",))
    ds.attrib["twist_model"] = s.twist_model
    defDim(ds, "x", s.nx)
    defDim(ds, "y", s.ny)
    defDim(ds, "t", s.N)
    ds.attrib["title"] = "4 phase decomposition"
    ds.attrib["casename"] = s.casename
    ds.attrib["twist"] = s.twist
    defVar(ds, "Lx", s.Lx, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "Ly", s.Ly, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "dx", s.dx, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "dy", s.dy, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "amplitude", s.amp, (), attrib=OrderedDict("units" => "m",))

    for var in ["eta", "phi"]
        v000 = ds_ep["$(var)000"] |> Array
        v090 = ds_ep["$(var)090"] |> Array
        v180 = ds_ep["$(var)180"] |> Array
        v270 = ds_ep["$(var)270"] |> Array
        v090h = ds_eph["$(var)090"] |> Array
        v270h = ds_eph["$(var)270"] |> Array
        for i in 1:4
            v = @match i begin
                1 => (v000 .- v090h .- v180 .+ v270h) ./ 4
                2 => (v000 .- v090 .+ v180 .- v270) ./ 4
                3 => (v000 .+ v090h .- v180 .- v270h) ./ 4
                4 => (v000 .+ v090 .+ v180 .+ v270) ./ 4
            end
            unit = @match var begin
                "eta" => "m"
                "phi" => "m²/s"
            end
            defVar(ds, "$var$i", v, ("t", "x", "y"), attrib=OrderedDict(
                "units" => unit,
            ))
        end
    end
    close(ds)
    @info @sprintf("Write complete for 4p.nc %s(%ddeg twist)", s.casename, s.twist)
end

function export_nc_kinematics(s; basedir=".")
    fname = joinpath(basedir, @sprintf("%s/%03ddeg/kinematics.nc", s.casename, s.twist))
    phases = [0, 90, 180, 270]
    if isfile(fname)
        @warn "NC file exists!"
        return
    end

    ds = NCDataset(fname, "c")
    ds.attrib["casename"] = s.casename
    ds.attrib["twist"] = s.twist
    ds.attrib["twist_model"] = s.twist_model
    defVar(ds, "Lx", s.Lx, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "Ly", s.Ly, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "dx", s.dx, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "dy", s.dy, (), attrib=OrderedDict("units" => "m",))
    defVar(ds, "dt", s.dt, (), attrib=OrderedDict("units" => "s",))
    defVar(ds, "amplitude", s.amp, (), attrib=OrderedDict("units" => "m",))

    for i in 1:4
        phase_str = @sprintf("%03d", phases[i])
        kf_path = joinpath(basedir, @sprintf("%s/%03ddeg/%s/1/Kinematics01.bin", s.casename, s.twist, phase_str))
        kf::KinematicsFile = open_Kinematics(kf_path)
        if i == 1
            defDim(ds, "x", kf.nx)
            defDim(ds, "y", kf.ny)
            defDim(ds, "z", kf.nz - 1)
            defDim(ds, "t", kf.nt)
        end
        defVar(ds, "eta$phase_str", kf.eta, ("t", "x", "y"), attrib=OrderedDict("units" => "m"))
        defVar(ds, "phi$phase_str", kf.phi, ("t", "x", "y"), attrib=OrderedDict("units" => "m²/s"))
        defVar(ds, "u$phase_str", kf.u[:, 2:end, :, :], ("t", "z", "x", "y"), attrib=OrderedDict("units" => "m/s"))
        defVar(ds, "v$phase_str", kf.v[:, 2:end, :, :], ("t", "z", "x", "y"), attrib=OrderedDict("units" => "m/s"))
        defVar(ds, "w$phase_str", kf.w[:, 2:end, :, :], ("t", "z", "x", "y"), attrib=OrderedDict("units" => "m/s"))
        defVar(ds, "uz$phase_str", kf.uz[:, 2:end, :, :], ("t", "z", "x", "y"), attrib=OrderedDict("units" => "m/s"))
        defVar(ds, "vz$phase_str", kf.vz[:, 2:end, :, :], ("t", "z", "x", "y"), attrib=OrderedDict("units" => "m/s"))
        defVar(ds, "wz$phase_str", kf.wz[:, 2:end, :, :], ("t", "z", "x", "y"), attrib=OrderedDict("units" => "m/s"))
    end

    close(ds)
    @info @sprintf("Write complete for kinematics.nc %s(%ddeg twist)", s.casename, s.twist)
end