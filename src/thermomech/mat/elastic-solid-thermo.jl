# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticSolidThermo

mutable struct ElasticSolidThermoIpState<:IpState
    env::ModelEnv
    σ::Array{Float64,1}
    ε::Array{Float64,1}
    V::Array{Float64,1}
    function ElasticSolidThermoIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.V = zeros(env.ndim)
        return this
    end
end


mutable struct ElasticSolidThermo<:Material
    E ::Float64 # Young's Modulus kN/m2
    nu::Float64
    k ::Float64 # thermal conductivity  w/m/k
    ρ ::Float64 # density Ton/m3
    cv::Float64 # Specific heat J/Ton/k
    α ::Float64 # thermal expansion coefficient  1/K or 1/°C


    function ElasticSolidThermo(prms::Dict{Symbol,Float64})
        return  ElasticSolidThermo(;prms...)
    end

    function ElasticSolidThermo(;E=1.0, nu=0.0, k=NaN, rho=NaN, cv=NaN, alpha=1.0)
        E<=0.0       && error("Invalid value for E: $E")
        !(0<=nu<0.5) && error("Invalid value for nu: $nu")
        isnan(k)     && error("Missing value for k")
        E<=0.0       && error("Invalid value for E: $E")
        cv<=0.0       && error("Invalid value for E: $E")
        0<=alpha<=1  || error("Invalid value for alpha: $alpha")
        this = new(E, nu, k, rho, cv, alpha)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::ElasticSolidThermo) = TMSolid

# Create a new instance of Ip data
new_ip_state(mat::ElasticSolidThermo, env::ModelEnv) = ElasticSolidThermoIpState(env)

function calcD(mat::ElasticSolidThermo, ipd::ElasticSolidThermoIpState)
    return calcDe(mat.E, mat.nu, ipd.env.modeltype) # function calcDe defined at elastic-solid.jl
end

function calcK(mat::ElasticSolidThermo, ipd::ElasticSolidThermoIpState) # Hydraulic conductivity matrix
    if ipd.env.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end

function stress_update(mat::ElasticSolidThermo, ipd::ElasticSolidThermoIpState, Δε::Array{Float64,1}, Δuθ::Float64, G::Array{Float64,1})
    De = calcD(mat, ipd)
    Δσ = De*Δε
    ipd.ε  += Δε
    ipd.σ  += Δσ
    K = calcK(mat, ipd)
    ipd.V   = -K*G
    return Δσ, ipd.V
end

function ip_state_vals(mat::ElasticSolidThermo, ipd::ElasticSolidThermoIpState)
    D = stress_strain_dict(ipd.σ, ipd.ε, ipd.env.ndim)

    D[:vx] = ipd.V[1]
    D[:vy] = ipd.V[2]
    if ipd.env.ndim==3
        D[:vz] = ipd.V[3]
    end

    return D
end
