# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export DruckerPrager

mutable struct DruckerPragerIpState<:IpState
    shared_data::SharedAnalysisData
    σ::Tensor2
    ε::Tensor2
    εpa::Float64
    Δγ::Float64
    function DruckerPragerIpState(shared_data::SharedAnalysisData=SharedAnalysisData()) 
        this = new(shared_data)
        this.σ   = zeros(6)
        this.ε   = zeros(6)
        this.εpa = 0.0
        this.Δγ  = 0.0
        this
    end
end

mutable struct DruckerPrager<:Material
    E::Float64
    ν::Float64
    α::Float64
    κ::Float64
    H::Float64
    ρ::Float64

    function DruckerPrager(prms::Dict{Symbol,Float64})
        return DruckerPrager(;prms...)
    end

    function DruckerPrager(;E=NaN, nu=0.0, alpha=0.0, kappa=0.0, H=0.0, rho=0.0)
        @assert E>0.0
        @assert 0.0<=nu<0.5
        @assert alpha>=0.0
        @assert kappa>0.0
        @assert H>=0.0
        @assert rho>=0.0
        
        this    = new(E, nu, alpha, kappa, H, rho)
        return this 
    end
end

matching_elem_type(::DruckerPrager) = MechSolid

# Create a new instance of Ip data
new_ip_state(mat::DruckerPrager, shared_data::SharedAnalysisData) = DruckerPragerIpState(shared_data)

function nlE(fc::Float64, εc::Float64, ε::Array{Float64,1})
    εv = abs(sum(ε[1:3]))
    return 2*fc*(εc-εv)/εc^2
end

function set_state(ipd::DruckerPragerIpState; sig=zeros(0), eps=zeros(0))
    if length(sig)==6
        ipd.σ .= sig.*V2M
    else
        if length(sig)!=0; error("DruckerPrager: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε .= eps.*V2M
    else
        if length(eps)!=0; error("DruckerPrager: Wrong size for strain array: $eps") end
    end
end

function yield_func(mat::DruckerPrager, ipd::DruckerPragerIpState, σ::Tensor2)
    j1  = J1(σ)
    j2d = J2D(σ)
    α,κ = mat.α, mat.κ
    H   = mat.H
    εpa = ipd.εpa
    return α*j1 + √j2d - κ - H*εpa
end

function calcD(mat::DruckerPrager, ipd::DruckerPragerIpState)
    α   = mat.α
    H   = mat.H
    De  = calcDe(mat.E, mat.ν, ipd.shared_data.model_type)

    if ipd.Δγ==0.0
        return De
    end

    j2d = J2D(ipd.σ)
    if j2d != 0.0
        s  = dev(ipd.σ) 
        su = s/norm(s)
        V  = α*tI + su/√2 # df/dσ
        N  = V
        Nu = N/norm(N)
    else # apex
        Nu = 1.0/√3.0*tI
        V  = Nu
    end

    return De - inner(De,Nu) ⊗ inner(V,De) / (inner(V,De,Nu) + H)
end

function stress_update(mat::DruckerPrager, ipd::DruckerPragerIpState, Δε::Tensor2)
    σini = ipd.σ
    De   = calcDe(mat.E, mat.ν, ipd.shared_data.model_type)
    σtr  = ipd.σ + inner(De, Δε)
    ftr  = yield_func(mat, ipd, σtr)

    if ftr < 1.e-8
        # elastic
        ipd.Δγ = 0.0
        ipd.σ  = σtr
    else
        # plastic 
        K, G  = mat.E/(3.0*(1.0-2.0*mat.ν)), mat.E/(2.0*(1.0+mat.ν))
        α, H  = mat.α, mat.H
        n     = 1.0/√(3.0*α*α+0.5)
        j1tr  = J1(σtr)
        j2dtr = J2D(σtr)

        if √j2dtr - ipd.Δγ*n*G > 0.0 # conventional return
            ipd.Δγ = ftr/(9*α*α*n*K + n*G + H)
            j1     = j1tr - 9*ipd.Δγ*α*n*K
            m      = 1.0 - ipd.Δγ*n*G/√j2dtr
            ipd.σ  = m*dev(σtr) + j1/3.0*tI
        else # return to apex
            κ      = mat.κ
            ipd.Δγ = (α*j1tr-κ-H*ipd.εpa)/(3*√3*α*K + H)
            j1     = j1tr - 3*√3*ipd.Δγ*K
            ipd.σ  = j1/3.0*tI
        end

        ipd.εpa += ipd.Δγ

    end

    ipd.ε += Δε
    Δσ     = ipd.σ - σini
    return Δσ
end


function ip_state_vals(mat::DruckerPrager, ipd::DruckerPragerIpState)
    ndim  = ipd.shared_data.ndim
    σ, ε  = ipd.σ, ipd.ε
    j1    = tr(σ)
    srj2d = √J2D(σ)

    D = stress_strain_dict(σ, ε, ndim)
    D[:epa]   = ipd.εpa
    D[:j1]    = j1
    D[:srj2d] = srj2d

    return D
end
