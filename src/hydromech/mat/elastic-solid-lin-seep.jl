# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticSolidLinSeep

mutable struct ElasticSolidLinSeepIpState<:IpState
    shared_data::SharedAnalysisData
    σ::Tensor2
    ε::Tensor2
    V::Array{Float64,1}
    uw::Float64
    function ElasticSolidLinSeepIpState(shared_data::SharedAnalysisData=SharedAnalysisData()) 
        this = new(shared_data)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.V = zeros(shared_data.ndim)
        this.uw = 0.0
        return this
    end
end


mutable struct ElasticSolidLinSeep<:Material
    E ::Float64
    nu::Float64
    k ::Float64
    gw::Float64

    function ElasticSolidLinSeep(prms::Dict{Symbol,Float64})
        return  ElasticSolidLinSeep(;prms...)
    end

    function ElasticSolidLinSeep(;E=1.0, nu=0.0, k=NaN, gw=NaN)
        E<=0.0       && error("Invalid value for E: $E")
        !(0<=nu<0.5) && error("Invalid value for nu: $nu")
        isnan(k)     && error("Missing value for k")
        isnan(gw)    && error("Missing value for gw")
        !(gw>0)      && error("Invalid value for gw: $gw")
        this    = new(E, nu, k, gw)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::ElasticSolidLinSeep) = HMSolid

# Create a new instance of Ip data
new_ip_state(mat::ElasticSolidLinSeep, shared_data::SharedAnalysisData) = ElasticSolidLinSeepIpState(shared_data)

function set_state(ipd::ElasticSolidLinSeepIpState; sig=zeros(0), eps=zeros(0))
    sq2 = √2.0
    mdl = [1, 1, 1, sq2, sq2, sq2]
    if length(sig)==6
        ipd.σ .= sig.*mdl
    else
        if length(sig)!=0; error("ElasticSolidLinSeep: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε .= eps.*mdl
    else
        if length(eps)!=0; error("ElasticSolidLinSeep: Wrong size for strain array: $eps") end
    end
end

function calcD(mat::ElasticSolidLinSeep, ipd::ElasticSolidLinSeepIpState)
    return calcDe(mat.E, mat.nu, ipd.shared_data.model_type) # function calcDe defined at elastic-solid.jl
end

function calcK(mat::ElasticSolidLinSeep, ipd::ElasticSolidLinSeepIpState) # Hydraulic conductivity matrix
    if ipd.shared_data.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end

function stress_update(mat::ElasticSolidLinSeep, ipd::ElasticSolidLinSeepIpState, Δε::Tensor2, Δuw::Float64, G::Array{Float64,1})
    De = calcD(mat, ipd)
    Δσ = De*Δε
    ipd.ε  += Δε
    ipd.σ  += Δσ
    K = calcK(mat, ipd)
    ipd.V   = -K*G
    ipd.uw += Δuw
    return Δσ, ipd.V
end

function ip_state_vals(mat::ElasticSolidLinSeep, ipd::ElasticSolidLinSeepIpState)
    D = stress_strain_dict(ipd.σ, ipd.ε, ipd.shared_data.ndim)

    D[:vx] = ipd.V[1]
    D[:vy] = ipd.V[2]
    if ipd.shared_data.ndim==3
        D[:vz] = ipd.V[3]
    end

    return D
end
