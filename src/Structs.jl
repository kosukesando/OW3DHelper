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