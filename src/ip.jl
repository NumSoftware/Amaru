# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

import Base.maximum
import Base.minimum
import Base.sort

"""
`IpState`

Abstract type for objects to store the state at integration points.
"""
abstract type IpState
    #env::ModelEnv
    #other data
end

function clone(src::IpState)
    dst = deepcopy(src)
    dst.env= src.env# keep original analyses data
    return dst
end

function Base.copy(src::IpState)
    T = typeof(src)
    dst = ccall(:jl_new_struct_uninit, Any, (Any,), T)
    names = fieldnames(T)
    for name in names
        val = getfield(src, name)
        if hasmethod(copy, (typeof(val),))
            setfield!(dst, name, copy(val))
        else
            setfield!(dst, name, val)
        end
    end
    return dst
end


function Base.copyto!(dst::IpState, src::IpState)
    names = fieldnames(typeof(src))
    for name in names
        val = getfield(src, name)
        if isbits(val)
            setfield!(dst, name, val)
        elseif typeof(val)<:AbstractArray
            copyto!(getfield(dst,name), val)
        else
            setfield!(dst, name, val)
            #error("copyto!(::IpState, ::IpState): unsupported field type: $(typeof(val))")
        end
    end
end


"""
`Ip(R, w)`

Creates an `Ip` object that represents an Integration Point in finite element analyses.
`R` is a vector with the integration point local coordinates and `w` is the corresponding integration weight.
"""
mutable struct Ip
    R    ::Array{Float64,1}
    w    ::Float64
    X    ::Array{Float64,1}
    id   ::Int
    tag  ::String
    owner::Any    # Element
    data ::IpState  # Ip current state

    function Ip(R::Array, w::Float64)
        this     = new(vec(R), w)
        this.X   = zeros(3)
        this.tag = ""
        this.owner = nothing
        return this
    end
end

# The functions below can be used in conjuntion with sort
get_x(ip::Ip) = ip.X[1]
get_y(ip::Ip) = ip.X[2]
get_z(ip::Ip) = ip.X[3]


"""
`ip_vals`

Returns a dictionary with keys and vals for the integration point `ip`.
"""
function ip_vals(ip::Ip)
    coords = Dict( :x => ip.X[1], :y => ip.X[2], :z => ip.X[3] )
    vals   = ip_state_vals(ip.owner.mat, ip.data)
    return merge(coords, vals)
end


# Index operator for a ip collection using expression
function Base.getindex(ips::Array{Ip,1}, filter_ex::Expr)
    R = Ip[]
    for ip in ips
        x, y, z = ip.X
        eval_arith_expr(filter_ex, x=x, y=y, z=z) && push!(R, ip)
    end
    return R
end


function Base.getindex(ips::Array{Ip,1}, s::String)
    return [ ip for ip in ips if ip.tag==s ]
end


# Get the maximum value of a given coordinate for the whole collection of ips
function maximum(ips::Array{Ip,1}, dir::Symbol)
    idx = findfisrt((:x, :y, :z), dir)
    maximum([ip.X[idx] for ip in ips])
end

function minimum(ips::Array{Ip,1}, dir::Symbol)
    idx = findfisrt((:x, :y, :z), dir)
    minimum([ip.X[idx] for ip in ips])
end

