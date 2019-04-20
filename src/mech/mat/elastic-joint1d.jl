# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticJoint1D

mutable struct Joint1DIpState<:IpState
    analysis_data::AnalysisData
    ndim::Int
    σ ::Array{Float64,1}
    u ::Array{Float64,1}
    function Joint1DIpState(analysis_data::AnalysisData=AnalysisData())
        this = new(analysis_data)
        ndim = analysis_data.ndim
        this.σ = zeros(ndim)
        this.u = zeros(ndim)
        return this
    end
end


mutable struct ElasticJoint1D<:Material
    ks::Float64
    kn::Float64
    h ::Float64    # section perimeter

    function ElasticJoint1D(prms::Dict{Symbol,Float64})
        return  ElasticJoint1D(;prms...)
    end

    function ElasticJoint1D(;ks=NaN, kn=NaN, h=NaN, A=NaN, dm=NaN)
        # A : section area
        # dm: section diameter
        # h : section perimeter
        ks>=0 || error("ks should be greater than zero")
        kn>=0 || error("kn should be greater than zero")
        (h>0 || A>0 || dm>0) || error("perimeter h, section area A or diameter dm should be provided")

        if isnan(h) 
            if A>0
                h = 2.0*(A*pi)^0.5
            else
                h = pi*dm
            end
        end
        @assert h>0

        this = new(ks, kn, h)
        return this
    end
end


# Returns the element type that works with this material
matching_elem_type(::ElasticJoint1D) = MechJoint1D

# Create a new instance of Ip data
new_ip_state(mat::ElasticJoint1D, analysis_data::AnalysisData) = Joint1DIpState(analysis_data)

function calcD(mat::ElasticJoint1D, ipd::Joint1DIpState)
    ks = mat.ks
    kn = mat.kn
    if ipd.analysis_data.ndim==2
        return [  ks  0.0 
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function stress_update(mat::ElasticJoint1D, ipd::Joint1DIpState, Δu)
    D = calcD(mat, ipd)
    Δσ = D*Δu

    ipd.u .+= Δu
    ipd.σ .+= Δσ
    return Δσ
end

function ip_state_vals(mat::ElasticJoint1D, ipd::Joint1DIpState)
    return OrderedDict(
      :ur   => ipd.u[1] ,
      :tau  => ipd.σ[1] )
end
