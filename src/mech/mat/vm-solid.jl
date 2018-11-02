export VonMises

mutable struct VonMisesIpState<:IpState
    shared_data::SharedAnalysisData
    σ::Tensor2
    ε::Tensor2
    εpa::Float64
    Δγ::Float64
    function VonMisesIpState(shared_data::SharedAnalysisData=SharedAnalysisData()) 
        this = new(shared_data)
        this.σ   = zeros(6)
        this.ε   = zeros(6)
        this.εpa = 0.0
        this.Δγ  = 0.0
        this
    end
end

mutable struct VonMises<:Material
    E ::Float64
    ν ::Float64
    σy::Float64
    H ::Float64
    #De::Tensor4

    function VonMises(prms::Dict{Symbol,Float64})
        return VonMises(;prms...)
    end

    function VonMises(;E=NaN, nu=0.0, σy=0.0, H=0.0)
        @assert E>0.0
        @assert 0.0<=nu<0.5
        @assert σy>0.0
        @assert H>=0.0

        this    = new(E, nu, σy, H)
        #this.De = calcDe(E, nu)
        return this 
    end
end

matching_elem_type(::VonMises) = MechSolid

# Create a new instance of Ip data
new_ip_state(mat::VonMises, shared_data::SharedAnalysisData) = VonMisesIpState(shared_data)

function set_state(ipd::VonMisesIpState; sig=zeros(0), eps=zeros(0))
    if length(sig)==6
        ipd.σ .= sig.*V2M
    else
        if length(sig)!=0; error("VonMises: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε .= eps.*V2M
    else
        if length(eps)!=0; error("VonMises: Wrong size for strain array: $eps") end
    end
end

function yield_func(mat::VonMises, ipd::VonMisesIpState, σ::Tensor2)
    j1  = J1(σ)
    j2d = J2D(σ)
    σy  = mat.σy 
    H   = mat.H
    εpa = ipd.εpa
    return √(3*j2d) - σy - H*εpa
end

function calcD(mat::VonMises, ipd::VonMisesIpState)
    σy = mat.σy 
    H  = mat.H
    #De = mat.De
    De  = calcDe(mat.E, mat.ν, ipd.shared_data.model_type)

    if ipd.Δγ==0.0
        return De
    end

    j2d = J2D(ipd.σ)
    @assert j2d>0
    s  = dev(ipd.σ) 
    su = s/norm(s)
    V  = √(3/2)*su # df/dσ
    Nu = su

    return De - inner(De,Nu) ⊗ inner(V,De) / (inner(V,De,Nu) + H)
end

function stress_update(mat::VonMises, ipd::VonMisesIpState, Δε::Tensor2)
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
        H     = mat.H
        j1tr  = J1(σtr)
        j2dtr = J2D(σtr)

        @assert √j2dtr - ipd.Δγ*√2*G > 0.0
        ipd.Δγ = ftr/(√6*G + H)
        j1     = j1tr
        m      = 1.0 - ipd.Δγ*√2*G/√j2dtr
        ipd.σ  = m*dev(σtr) + j1/3.0*tI

        ipd.εpa += ipd.Δγ

    end

    ipd.ε += Δε
    Δσ     = ipd.σ - σini
    return Δσ
end

function ip_state_vals(mat::VonMises, ipd::VonMisesIpState)
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
