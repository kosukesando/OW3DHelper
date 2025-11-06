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

struct KinematicsFile
    dt::Float64
    nz::Int
    nx::Int
    ny::Int
    nt::Int
    t::Vector{Float64}
    x::Array{Float64}
    y::Array{Float64}
    eta::Array{Float64}
    phi::Array{Float64}
    u::Array{Float64}
    uz::Array{Float64}
    v::Array{Float64}
    vz::Array{Float64}
    w::Array{Float64}
    wz::Array{Float64}
end

struct KinematicsFileFull
    xbeg::Int
    xend::Int
    xstride::Int
    ybeg::Int
    yend::Int
    ystride::Int
    tbeg::Int
    tend::Int
    tstride::Int
    dt::Float64
    nz::Int
    nx::Int
    ny::Int
    nt::Int
    sigma::Vector{Float64}
    t::Vector{Float64}
    x::Array{Float64}
    y::Array{Float64}
    h::Array{Float64}
    hx::Array{Float64}
    hy::Array{Float64}
    eta::Array{Float64}
    etax::Array{Float64}
    etay::Array{Float64}
    phi::Array{Float64}
    u::Array{Float64}
    uz::Array{Float64}
    v::Array{Float64}
    vz::Array{Float64}
    w::Array{Float64}
    wz::Array{Float64}
end

struct PostProcessSetting
    casename::String
    amp::Float64
    phase::Float64
    twist::Int
    twist_model::String
    nt::Int
    nx::Int
    ny::Int
    dt::Float64
    dx::Float64
    dy::Float64
    Lx::Float64
    Ly::Float64
end
