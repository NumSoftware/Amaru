# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

__precompile__() 

"""
**Amaru.jl**

Amaru module implements functions and types to perform finite element analyses.

**Important data types**

Node, Element, Domain, Dof, Ip, NodeBC, FaceBC

**Important functions** 

set_mat, set_bc, clear_bc, solve!, save

"""
module Amaru
using Printf, Statistics, LinearAlgebra, SparseArrays, DelimitedFiles, Arpack
using JSON, DataStructures, Reexport
import DataStructures.OrderedDict, DataStructures.OrderedSet

#Using non-registered (jet) package FemMesh
#try
    #eval(:(using FemMesh))
#catch err
    #if isa(err, ArgumentError)
        #contains(err.msg, "FemMesh not found") && pkg.clone("https://github.com/NumSoftware/FemMesh")
    #else
        #error("Amaru: Error loading FemMesh package.")
    #end
#end

@reexport using FemMesh
import FemMesh.save # to be extended
import FemMesh.update! # to be extended
import FemMesh.get_x, FemMesh.get_y, FemMesh.get_z # to be extended
import FemMesh.tag!


# Debug
mutable struct DebugFlags
    enabled::Bool
    makeplots::Bool
end
const debug = DebugFlags(false, true)

# eye function
eye(n::Int64) = Array{Float64}(I,n,n)

# Tools module
include("tools/constants.jl")
include("tools/math.jl")
include("tools/linalg.jl")
include("tools/expr.jl")
include("tools/table.jl")
include("tools/tensors.jl")
include("tools/show.jl")
include("tools/stopwatch.jl")

# generic exports
export max, min, sort, reset, getindex, sort, copy!, show

# Fem module
include("model-env.jl")
export ModelEnv

include("node.jl")
export Node, Dof, add_dof, nodes_dof_vals

include("ip.jl")
export Ip, ip_vals, maximum, minimum, sort

include("material.jl")
export Material, read_prms

include("element.jl")
export Element
export get_nodes, get_ips, getcoords, getvals, get_map, elems_ip_vals

include("facet.jl")
export Facet, Face, Edge

include("tag.jl")
export tag!

include("bc.jl")
#export NodeBC, FaceBC, EdgeBC, apply_bc
export NodeBC, FaceBC, EdgeBC, ElemBC, apply_bc

include("logger.jl")
export Logger, GroupLogger
export NodeLogger, IpLogger, FaceLogger, EdgeLogger, NodeGroupLogger, IpGroupLogger
export update_logger!

include("domain.jl")
export Domain, SubDomain, reset!, setloggers!, datafields


# Mechanical module
include("mech/include.jl")

# Hydromechanical module
include("hydromech/include.jl")

# show function for Amaru types
for datatype in (:Dof, :Node, :Ip, :IpState, :Element, :Material, :BC, :Facet, :AbstractLogger, :Domain)
    eval( quote
        function Base.show(io::IO, obj::$datatype)
            print_field_values(io, obj)
        end

        function Base.show(io::IO, array::Array{<:$datatype,1})
            print_array_values(io, array)
        end
    end )
end


end#module
