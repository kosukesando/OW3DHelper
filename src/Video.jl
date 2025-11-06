function plot_raw2d(s, var; basedir=".", xlim=(-1000, 1000), ylim=(-1000, 1000))
    @assert var in ["eta", "phi"]
    dir = joinpath(basedir, @sprintf("%s/%03ddeg", s.casename, s.twist))
    mktempdir() do tempdir
        Threads.@threads for i in 1:s.N
            local fname = joinpath(tempdir, @sprintf("%03d.png", i))
            local v = NCDataset(joinpath(dir, "ep.nc"))[@sprintf("%s%03d", var, s.phase)][i, :, :]
            local f = Figure(size=(11cm, 8cm), fontsize=12pt)
            local ax = Axis(
                f[1, 1],
                title="Time=$(i)s",
                xlabel="x(m)",
                ylabel="y(m)",
            )
            colsize!(f.layout, 1, Aspect(1, 1.0))
            xlims!(ax, xlim)
            ylims!(ax, ylim)
            local im = image!(
                ax,
                (-s.Lx / 2, s.Lx / 2),
                (-s.Ly / 2, s.Ly / 2),
                v,
                colorrange=(-s.amp * 1.2, s.amp * 1.2),
                colormap=:seismic,
                highclip=:black,
            )
            Colorbar(f[1, 2], im)
            save(fname, f)
        end
        output = joinpath(dir, @sprintf("videos/%s_2d_%03d.mp4", var, s.phase))
        export_video(tempdir, output)
    end
end

function plot_raw3d(s::PostProcessSetting, var; basedir=".", xlim=(-1000, 1000), ylim=(-1000, 1000))
    @assert var in ["eta", "phi"]
    dir = joinpath(basedir, @sprintf("%s/%03ddeg", s.casename, s.twist))
    mktempdir() do tempdir
        Threads.@threads for i in 1:s.nt
            local fname = joinpath(tempdir, @sprintf("%03d.png", i))
            local v = NCDataset(joinpath(dir, "ep.nc"))[@sprintf("%s%03d", var, s.phase)][i, :, :]
            local f = Figure(size=(11cm, 8cm), fontsize=12pt)
            local ax = Axis3(f[1, 1],
                title="Time=$(i)s",
                xlabel="x (m)",
                ylabel="y (m)",
                zlabel="z (m)",
                azimuth=-pi / 4,
                elevation=pi / 3,
            )
            xlims!(ax, -1000, 1000)
            ylims!(ax, -1000, 1000)
            zlims!(ax, -s.amp, s.amp)
            local sf = surface!(ax,
                -s.Lx/2:s.dx:s.Lx/2,
                -s.Ly/2:s.dy:s.Ly/2,
                v,
                colormap=:deepsea,
                colorrange=(-s.amp * 1.2, s.amp * 1.2),
            )
            Colorbar(f[1, 2], sf)
            save(fname, f)
        end
        output = joinpath(dir, @sprintf("videos/%s_3d_%03d.mp4", var, s.phase))
        export_video(tempdir, output)
    end
end
function plot_4p(s, var, combination; basedir=".")
    @assert var in ["eta", "phi"]
    @assert combination in 1:4
    if combination == 1
        mult = 1
    else
        mult = 0.1
    end
    vlim = s.amp * mult
    dir = @sprintf("%s/%03ddeg", s.casename, s.twist)
    mktempdir() do tempdir
        Threads.@threads for i in 1:s.N
            local fname = joinpath(tempdir, @sprintf("%03d.png", i))
            local fp = NCDataset(joinpath(dir, "4p.nc"))[@sprintf("%s%d", var, combination)][i, :, :]
            local f = Figure(size=(11cm, 8cm), fontsize=12pt)
            local ax = Axis3(f[1, 1],
                title="Time=$(i)s",
                xlabel="x (m)",
                ylabel="y (m)",
                zlabel="z (m)",
                azimuth=-pi / 4,
                elevation=pi / 3,
            )
            xlims!(ax, -1000, 1000)
            ylims!(ax, -1000, 1000)
            zlims!(ax, -vlim, vlim)
            local sf = surface!(ax,
                -s.Lx/2:s.dx:s.Lx/2,
                -s.Ly/2:s.dy:s.Ly/2,
                fp,
                colormap=:deepsea,
                colorrange=(-vlim, vlim),
            )
            Colorbar(f[1, 2], sf)
            save(fname, f)
        end
        output = joinpath(dir, "videos/4p$combination.mp4")
        export_video(tempdir, output)
    end
end

function export_video(tempdir, output)
    mkpath(dirname(output))
    command = `ffmpeg -i $(tempdir)/%03d.png -vcodec libx264 -b:v 1600k -pix_fmt yuv420p -y $(output)`
    try
        run(command)
    catch e
        cp(tempdir, splitext(output)[1], force=true)
    end
end
