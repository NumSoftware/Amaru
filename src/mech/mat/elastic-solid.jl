# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticSolid

mutable struct ElasticSolidIpState<:IpState
    shared_data::SharedAnalysisData
    #σ::Array{Float64,1}
    #ε::Array{Float64,1}
    σ::Tensor2
    ε::Tensor2
    function ElasticSolidIpState(shared_data::SharedAnalysisData=SharedAnalysisData()) 
        this = new(shared_data)
        this.σ = zeros(6)
        this.ε = zeros(6)
        return this
    end
end


mutable struct ElasticSolid<:Material
    E ::Float64
    nu::Float64
    ρ::Float64

    function ElasticSolid(prms::Dict{Symbol,Float64})
        return  ElasticSolid(;prms...)
    end

    function ElasticSolid(;E=1.0, nu=0.0, rho=0.0)
        E<=0.0       && error("Invalid value for E: $E")
        !(0<=nu<0.5) && error("Invalid value for nu: $nu")
        rho<0.0      && error("Invalid value for rho: $rho")
        this = new(E, nu, rho)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::ElasticSolid) = MechSolid

# Create a new instance of Ip data
new_ip_state(mat::ElasticSolid, shared_data::SharedAnalysisData) = ElasticSolidIpState(shared_data)

function set_state(ipd::ElasticSolidIpState; sig=zeros(0), eps=zeros(0))
    sq2 = √2.0
    mdl = [1, 1, 1, sq2, sq2, sq2]
    if length(sig)==6
        ipd.σ .= sig.*mdl
    else
        if length(sig)!=0; error("ElasticSolid: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε .= eps.*mdl
    else
        if length(eps)!=0; error("ElasticSolid: Wrong size for strain array: $eps") end
    end
end

function calcDe(E::Number, ν::Number, model_type::Symbol)
    if model_type==:plane_stress
        c = E/(1.0-ν^2)
        return Tensor4( [
            c    c*ν   0.0  0.0  0.0  0.0
            c*ν  c     0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0 
            0.0  0.0   0.0  0.0  0.0  c*(1.0-ν) ] )
        # ezz = -ν/E*(sxx+syy)
    else
        c = E/((1.0+ν)*(1.0-2.0*ν))
        return Tensor4( [
            c*(1-ν) c*ν     c*ν     0.0         0.0         0.0
            c*ν     c*(1-ν) c*ν     0.0         0.0         0.0
            c*ν     c*ν     c*(1-ν) 0.0         0.0         0.0
            0.0     0.0     0.0     c*(1-2*ν)   0.0         0.0
            0.0     0.0     0.0     0.0         c*(1-2*ν)   0.0
            0.0     0.0     0.0     0.0         0.0         c*(1-2*ν) ] )
    end
end

function calcD(mat::ElasticSolid, ipd::ElasticSolidIpState)
    return calcDe(mat.E, mat.nu, ipd.shared_data.model_type)
end

function stress_update(mat::ElasticSolid, ipd::ElasticSolidIpState, dε::Tensor2)
    De = calcDe(mat.E, mat.nu, ipd.shared_data.model_type)
    dσ = De*dε
    ipd.ε += dε
    ipd.σ += dσ
    return dσ
end

function ip_state_vals(mat::ElasticSolid, ipd::ElasticSolidIpState)
    return stress_strain_dict(ipd.σ, ipd.ε, ipd.shared_data.ndim)
end
