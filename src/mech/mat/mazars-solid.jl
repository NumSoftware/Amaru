# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export Mazars

mutable struct MazarsIpState<:IpState
    shared_data::SharedAnalysisData
    σ::Tensor2
    ε::Tensor2
    φt::Float64
    φc::Float64
    φ::Float64  # damage
    ε̅max::Float64
    function MazarsIpState(shared_data::SharedAnalysisData=SharedAnalysisData()) 
        this = new(shared_data)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.φ = 0.0
        this.φc = 0.0
        this.φt = 0.0
        this.ε̅max = 0.0
        return this
    end
end

mutable struct Mazars<:Material
    E ::Float64
    nu::Float64
    ε̅0::Float64
    At::Float64
    Bt::Float64
    Ac::Float64
    Bc::Float64
    ρ::Float64
    De::Tensor4
    invDe::Tensor4

    function Mazars(prms::Dict{Symbol,Float64})
        return Mazars(;prms...)
    end

    function Mazars(;E=NaN, nu=0.0, eps0=NaN, At=NaN, Bt=NaN, Ac=NaN, Bc=NaN, rho=0.0)
        @assert E>0.0
        @assert 0.0<=nu<0.5
        @assert At>0.0
        @assert Ac>0.0
        @assert Bt>0.0
        @assert Bc>0.0
        @assert rho>=0.0
        @assert eps0>0

        this     = new(E, nu, eps0, At, Bt, Ac, Bc, rho)
        this.De  = calcDe(E, nu, :general)
        this.invDe  = inv(this.De)
        return this 
    end
end

# Returns the element type that works with this material model
matching_elem_type(::Mazars) = MechSolid

# Create a new instance of Ip data
new_ip_state(mat::Mazars, shared_data::SharedAnalysisData) = MazarsIpState(shared_data)

function set_state(ipd::MazarsIpState; sig=zeros(0), eps=zeros(0))
    if length(sig)==6
        ipd.σ .= sig.*V2M
    else
        if length(sig)!=0; error("Mazars: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε .= eps.*V2M
    else
        if length(eps)!=0; error("Mazars: Wrong size for strain array: $eps") end
    end
end


# special functions
pos(x) = (abs(x)+x)/2.0
neg(x) = (-abs(x)+x)/2.0


function calcD(mat::Mazars, ipd::MazarsIpState)
    # There is something wrong with the derivatives here

    # Equivalent strain scalar
    εp = eigenvals(ipd.ε)
    ε̅ = norm(pos.(εp))
    ε̅ == 0.0 && (ε̅ += 1e-15)
    ε̅max = max(ipd.ε̅max, mat.ε̅0)

    if ε̅<ε̅max
        #@show "elastic"
        return (1.0 - ipd.φ)*mat.De
    else
        #@show "plastic"
        #return (1.0 - ipd.φ)*mat.De
        return mat.De

        # Principal stresses and principal directions
        σp, V = eigen(ipd.σ)
        σp = [ σp; zeros(3) ]
        # Eigen vectors
        p1 = V[:, 1]
        p2 = V[:, 2]
        p3 = V[:, 3]
        P  = [ matrix2Mandel(p1*p1'), matrix2Mandel(p2*p2'), matrix2Mandel(p3*p3') ]

        # Tensile and compression tensors
        σt = pos.(σp)
        σc = neg.(σp)
        Dei= inv(mat.De)
        εt = Dei*σt
        εc = Dei*σc
        εv = sum(pos.(εt)) + sum(neg.(εc))

        εv +=  1e-15 # avoid division by zero

        # Tensile and compression damage weights
        αt = clamp(sum(pos.(εt))/εv, 0.0, 1.0)
        αc = clamp(sum(neg.(εc))/εv, 0.0, 1.0)

        # Constitutive matrix calculation
        dφtdε̅ = (1.0-mat.At)*mat.ε̅0*ε̅^-2 + mat.At*mat.Bt*exp( -mat.Bt*(ε̅-mat.ε̅0) )
        dφcdε̅ = (1.0-mat.Ac)*mat.ε̅0*ε̅^-2 + mat.Ac*mat.Bc*exp( -mat.Bc*(ε̅-mat.ε̅0) )
        dε̅dε  = sum( pos(εp[i])*(1+sign(εp[i]))*P[i] for i=1:3 ) / (2*ε̅)
        dφdε  = (αt*dφtdε̅ + αc*dφcdε̅) * dε̅dε

        #@show dφdε'*mat.De
        #D     = (1.0 - ipd.φ)*mat.De - (dφdε'*mat.De)'*ipd.ε'
        D     = (1.0 - ipd.φ)*mat.De - dφdε*(mat.De*ipd.ε)'
        return D
    end
end


function stress_update(mat::Mazars, ipd::MazarsIpState, Δε::Tensor2)
    σini  = ipd.σ
    ipd.ε = ipd.ε + Δε

    E  = mat.E
    nu = mat.nu

    # Principal stresses tensor
    εp = eigenvals(ipd.ε)

    # Equivalent strain scalar
    ε̅ = norm(pos.(εp))
    ε̅ == 0.0 && (ε̅ += 1e-15)
    ipd.ε̅max = max(ipd.ε̅max, mat.ε̅0)

    if ε̅ < ipd.ε̅max  # linear-elastic increment
        #@show "Elast"
        ipd.σ = (1.0 - ipd.φ)*mat.De*ipd.ε
    else # increment with damage: A previous elastic step is desired
        #@show "Plast"
        ipd.ε̅max = ε̅

        # Principal stresses and principal directions
        σp = eigenvals(ipd.σ)
        σp = [ σp; zeros(3) ]

        # Damage calculation
        ipd.φt = 1.0 - (1-mat.At)*mat.ε̅0/ε̅ - mat.At/exp(mat.Bt*(ε̅-mat.ε̅0))
        ipd.φc = 1.0 - (1-mat.Ac)*mat.ε̅0/ε̅ - mat.Ac/exp(mat.Bc*(ε̅-mat.ε̅0))

        # Tensile and compression tensors
        σt = pos.(σp)
        σc = neg.(σp)
        εt = mat.invDe*σt
        εc = mat.invDe*σc

        #εv = sum(pos(εt[i]) + pos(εc[i]) for i=1:3)
        εv = sum(pos.(εt)) + sum(neg.(εc))
        εv == 0.0 && (εv += 1e-15)

        # Tensile and compressive damage weights
        αt = clamp(sum(pos.(εt))/εv, 0.0, 1.0)
        αc = clamp(sum(neg.(εc))/εv, 0.0, 1.0)

        # Damage variable
        φ = αt*ipd.φt + αc*ipd.φc
        ipd.φ = clamp(φ, ipd.φ, 0.999)

        # Total stress and stress increment
        ipd.σ = (1.0 - ipd.φ)*mat.De*ipd.ε
    end

    Δσ    = ipd.σ - σini
    return Δσ 
end

function ip_state_vals(mat::Mazars, ipd::MazarsIpState)
    ndim  = ipd.shared_data.ndim
    σ, ε  = ipd.σ, ipd.ε

    D = stress_strain_dict(σ, ε, ndim)
    D[:dam]  = ipd.φ
    D[:damt] = ipd.φt
    D[:damc] = ipd.φc
    D[:eq]   = ipd.ε̅max

    return D
end
