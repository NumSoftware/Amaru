# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export Orthotropic


mutable struct FixedPlane
    V::Array{Float64,1} # direction
    failed::Bool
    active::Bool # active failure
    εf::Float64  # strain at failure
end


mutable struct OrthotropicIpState<:IpState
    shared_data::SharedAnalysisData
    σ::Tensor2
    ε::Tensor2
    fmax::Float64

    unloading::Bool  # flag for loading/unloading conditions
    crushed::Bool    # flag for crushing in compression

    fplanes::Array{FixedPlane,1}
    ep3max::Float64
    sp3max::Float64

    function OrthotropicIpState(shared_data::SharedAnalysisData=SharedAnalysisData()) 
        this = new(shared_data)
        this.σ  = zeros(6)
        this.ε  = zeros(6)
        this.fmax = 0.0
        this.unloading = true
        this.crushed   = false
        this.fplanes = [ ]
        this.ep3max = 0.0
        this.sp3max = 0.0
        this
    end
end

mutable struct Orthotropic<:Material
    E0::Float64  # Young's modulus
    ν ::Float64  # Poisson ratio
    ft::Float64  # compression strenght
    fc::Float64  # tensile strenght
    fu::Float64  # ultimate compression strenght
    εc::Float64  # strain corresponding to compression strenght
    εu::Float64  # strain corresponding to ultimate compression strenght
    α ::Float64  # hydrostatic multiplier for the loading function (similar do DP model)
    ηn::Float64  # reduction factor for normal stiffness components
    ηs::Float64  # reduction factor for shear stiffness components

    function Orthotropic(prms::Dict{Symbol,Float64})
        return Orthotropic(;prms...)
    end

    function Orthotropic(;E=NaN, nu=0.0, ft=0.0, fc=0.0, fu=0.0, epsc=0.0, epsu=0.0, alpha=0.4, etan=0.001, etas=0.5)
        @assert E>0.0       
        @assert 0.0<=nu<0.5 
        @assert fc>0.0      
        @assert ft>0.0      
        @assert fu<fc       
        @assert alpha>0     
        @assert epsc>0      
        @assert epsu>epsc  # εu>εc
        @assert etan>0     # ηn
        @assert etas>0     # ηs
        
        this = new(E, nu, ft, -fc, -fu, -epsc, -epsu, alpha, etan, etas)
        return this 
    end
end

matching_elem_type(::Orthotropic) = MechSolid

# Create a new instance of Ip data
new_ip_state(mat::Orthotropic, shared_data::SharedAnalysisData) = OrthotropicIpState(shared_data)

function set_state(ipd::OrthotropicIpState; sig=zeros(0), eps=zeros(0))
    if length(sig)==6
        ipd.σ[:] = sig.*V2M
    else
        if length(sig)!=0; error("Orthotropic: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε[:] = eps.*V2M
    else
        if length(eps)!=0; error("Orthotropic: Wrong size for strain array: $eps") end
    end
end


function loading_func(mat::Orthotropic, σ::Tensor2)
    j2d = J2D(σ)
    j1  = J1(σ)
    return √j2d
    #return √j2d + mat.α*j1
end


function sigma(mat::Orthotropic, εp::Float64)
    # εp : principal strain
    # σp : principal stress

    E0 = mat.E0
    Eu = mat.fu/mat.εu
    Es = mat.fc/mat.εc
    p  = mat.εu/mat.εc
    A = ( E0/Eu + (p^3-2*p^2)*E0/Es - (2*p^3-3*p^2+1) ) / (p^3-2*p^2+p)
    B = ( 2*E0/Es - 3 ) - 2*A
    C = ( 2 - E0/Es ) + A
    ξ = εp/mat.εc
    return (E0/Es*ξ) / ( 1 + A*ξ + B*ξ^2 + C*ξ^3) * mat.fc
end


function uniaxial_young_modulus(mat::Orthotropic, εp::Float64, γ1::Float64=1.0)
    # εp : principal strain
    # σp : principal stress

    gamma = 1.0 # factor for strains
    εp = γ1*gamma*εp
    εu = γ1*gamma*mat.εu

    E0 = mat.E0
    Eu = mat.fu/mat.εu
    Es = mat.fc/mat.εc
    p  = mat.εu/mat.εc
    A = ( E0/Eu + (p^3-2*p^2)*E0/Es - (2*p^3-3*p^2+1) ) / (p^3-2*p^2+p)
    B = 2*E0/Es - 3 - 2*A
    C = 2 - E0/Es + A
    ξ = εp/mat.εc
    E = E0*(1 - B*ξ^2 - 2*C*ξ^3) / ( 1 + A*ξ + B*ξ^2 + C*ξ^3)^2
    return E
end

function gamma1(mat::Orthotropic, σp1::Float64, σp2::Float64)
    # Calculate the amplifying factor γ1 for compression strenght due to confinning stresses σp1 and σp2

    λ = σp1/mat.fc
    β = σp2/mat.fc

    λ<0 || β<0 && return 1.0

    p1 = [ 1.2, 1.2 ]
    p2 = [ 0.8, 1.25 ]
    p3 = [ 0.0, 1.0 ]

    p1λ = p1 .+ 0.55*λ/0.25*normalize([1.0 , 1.0])
    p2λ = p2 .+ 0.50*λ/0.25*normalize([0.75, 1.0])
    p3λ = p3 .+ 0.50*λ/0.25*normalize([0.75, 1.0])

    gamma1 = p1λ[2] * (β-p2λ[1]) * (β-p3λ[1]) / ( (p1λ[1]-p2λ[1]) * (p1λ[1]-p3λ[1]) ) +
             p2λ[2] * (β-p1λ[1]) * (β-p3λ[1]) / ( (p2λ[1]-p1λ[1]) * (p2λ[1]-p3λ[1]) ) +
             p3λ[2] * (β-p1λ[1]) * (β-p2λ[1]) / ( (p3λ[1]-p1λ[1]) * (p3λ[1]-p2λ[1]) )

    gamma1 = max(gamma1, 1.0)

    return gamma1
end

function eigen_with_fixed_dir(σ::Tensor2, X::Array{Float64,1})
    # Finds two eigenvectors normal to X. Vector X is returned as the first direction
    
    # σ: Stress tensor
    # X: Fixed direction

    # find an arbitrary system aligned with X
    Q = [1., 0 , 0] # auxiliary vector
    if X==Q
        Q = [0., 1, 0 ]
    end
    Y = normalize(cross(X, Q))
    Z = normalize(cross(X, Y))
    W = [X Y Z] # arbitrary system

    # fourth order rotation tensor
    R = zeros(6,6)
    tensor_rot!(W, R)

    # new temporary tensor
    σt = R*σ
    σx = σt[1] # first "eigenvalue"

    lt, Vt = eigen( [ σt[2] σt[4]; σt[4] σt[3] ]  )

    # 2D to 3D
    Dt = eye(3)
    Dt[2:3,2:3] .= Vt # new matrix with orthogonal directions
    Dt[:,3] .= cross(Dt[:,1], Dt[:,2])
    #@show Dt[:,1]

    # new tensor
    V = W*Dt  # new directions in the xyz system
    L  = [ σx; lt ]  # new orthotropic stresses
    return L, V

end

function fixD(mat, D::Array{Float64,2}, active_fails::Array{Int,1})
    ηn = mat.ηn
    ηs = mat.ηs
    for i in active_fails
        D[i,1:3] *= ηn
        D[1:3,i] *= ηn
        if i==1
            D[5,5] *= ηs
            D[6,6] *= ηs
        elseif i==2
            D[4,4] *= ηs
            D[6,6] *= ηs
        elseif i==3
            D[4,4] *= ηs
            D[5,5] *= ηs
        end
    end

    return D
end

function find_dir(V, Vi)
    tol = 1e-5
    norm( V[:,1] - Vi) < tol && return 1
    norm( V[:,2] - Vi) < tol && return 2
    norm( V[:,3] - Vi) < tol && return 3
error("Plano não encontrado")
end

function principal_isotropic_moduli(mat::Orthotropic, εp::Array{Float64,1}, γ1::Float64)
    Ep1 = uniaxial_young_modulus(mat, εp[1], γ1)
    Ep2 = uniaxial_young_modulus(mat, εp[2], γ1)
    Ep3 = uniaxial_young_modulus(mat, εp[3], γ1)
    return Ep1, Ep2, Ep3
end

function principal_orthotropic_moduli(mat::Orthotropic, εp::Array{Float64,1}, Δεp::Array{Float64,1}, γ1::Float64)
    Ep1 = ( sigma(mat, εp[1]+Δεp[1]) - sigma(mat, εp[1]) ) / Δεp[1]
    Ep2 = ( sigma(mat, εp[2]+Δεp[2]) - sigma(mat, εp[2]) ) / Δεp[2]
    Ep3 = ( sigma(mat, εp[3]+Δεp[3]) - sigma(mat, εp[3]) ) / Δεp[3]
    Δεp[1] == 0 && (Ep1 = uniaxial_young_modulus(mat, εp[1], γ1))
    Δεp[2] == 0 && (Ep2 = uniaxial_young_modulus(mat, εp[2], γ1))
    Δεp[3] == 0 && (Ep3 = uniaxial_young_modulus(mat, εp[3], γ1))
    return Ep1, Ep2, Ep3
end

function orthotropic_moduli(Ep1::Float64, Ep2::Float64, Ep3::Float64, σp::Array{Float64,1})
    E12 = ( abs(σp[1])*Ep1 + abs(σp[2])*Ep2 ) / ( abs(σp[1]) + abs(σp[2]) )
    E23 = ( abs(σp[2])*Ep2 + abs(σp[3])*Ep3 ) / ( abs(σp[2]) + abs(σp[3]) )
    E13 = ( abs(σp[1])*Ep1 + abs(σp[3])*Ep3 ) / ( abs(σp[1]) + abs(σp[3]) )
    σp[1] + σp[2] == 0.0 && (E12=Ep1)
    σp[2] + σp[3] == 0.0 && (E23=Ep2)
    σp[1] + σp[3] == 0.0 && (E13=Ep3)
    return E12, E23, E13
end


function calcD(mat::Orthotropic, ipd::OrthotropicIpState)
    println("------------------------------------------")
    println("Caldulo de D")

    E = mat.E0
    ν = mat.ν

    nactive_fails = length(ipd.fplanes)

    nfplanes = length(ipd.fplanes)
    if nfplanes>0
        nactive_fails = sum( p.active for p in ipd.fplanes )
    else
        nactive_fails = 0
    end

    if !ipd.crushed && nactive_fails==0 # no fails
        println("Não está crushed e não tem planos de falha")
        if ipd.unloading # elastic regime
            println("Descarregando")
            println("De")
            D = calcDe(mat.E0, mat.ν, ipd.shared_data.model_type)
            return D
        else # loading
            println("Carregando")
            σp, V = eigen(ipd.σ)
            p  = sortperm(σp, rev=true)
            σp = σp[p] # ordered stresses
            V  = V[:,p]

            R = zeros(6,6)
            tensor_rot!(V, R)
            εp = R*ipd.ε # strain associated with σp

            γ1  = gamma1(mat, σp[1], σp[2])
            fcm = γ1*mat.fc

            # Use secant Young's moduli
            Ep1, Ep2, Ep3 = principal_isotropic_moduli(mat, εp, γ1)

            κ = 0.4
            if σp[3] >= κ*fcm # σc'  low compression
                println("Baixa compressão ou tração")
                println("D(isotrop)")
                if σp[1]!=0 || σp[2]!=0 || σp[3]!=0 
                    Et = ( abs(σp[1])*Ep1 + abs(σp[2])*Ep2 + abs(σp[3])*Ep3 ) / (abs(σp[1]) + abs(σp[2]) + abs(σp[3]))
                else
                    Et = mat.E0
                end
                D = calcDe(Et, mat.ν, ipd.shared_data.model_type)
                return D

            else # σp[3]<κ*fcmax high compression
                println("Alta compressão")
                println("D(ortho)")
                E12, E23, E13 = orthotropic_moduli(Ep1, Ep2, Ep3, σp)

                # Orthotropic D matrix
                # Notice that Amaru considers a general shear stress components (e.g. εxy)
                # and not the engineering definitions (e.g, γxy).
                Dp = 1/((1+ν)*(1-2*ν))*[ (1-ν)*Ep1   ν*E12       ν*E13      0.0          0.0          0.0 
                                         ν*E12       (1-ν)*Ep2   ν*E23      0.0          0.0          0.0 
                                         ν*E13       ν*E23      (1-ν)*Ep3   0.0          0.0          0.0 
                                         0.0         0.0         0.0        (1-2*ν)*E12  0.0          0.0
                                         0.0         0.0         0.0        0.0          (1-2*ν)*E13  0.0  
                                         0.0         0.0         0.0        0.0          0.0          (1-2*ν)*E23 ]
                D = R'*Dp*R # rotates tensor Dp to xyz system
                return D
            end
        end

    elseif ipd.crushed
        println("Crushed")
        println("D=ηn")
        D = eye(6)*mat.ηn
        return D
    else # nactive_fails>0
        println("Tração")
        if nfplanes==0
            println("Planos de falha = 0")
            σp, V = eigen(ipd.σ)
        elseif nfplanes==1
            println("Plano de falha definido = 1")
            σp, V = eigen_with_fixed_dir(ipd.σ, ipd.fplanes[1].V) # V1 should be the first column of V
        else # nfplanes == 3
            println("Planos de falha definido = 2 ou 3")
            V = [ ipd.fplanes[1].V ipd.fplanes[2].V ipd.fplanes[3].V ]
            R = zeros(6,6)
            tensor_rot!(V, R)
            σp = R*ipd.σ # stresses associated with V1, V2 and V3
            σp = σp[1:3]
        end

        # sort principal stresses
        p  = sortperm(σp, rev=true)
        σp = σp[p] # ordered stresses
        V  = V[:,p]
        R  = zeros(6,6)
        tensor_rot!(V, R)

        if ipd.unloading
        println("Descarregando")
            Dp = calcDe(mat.E0, mat.ν, ipd.shared_data.model_type)

        else # loading
            println("Carregando")
            εp  = R*ipd.ε # strain associated with σp
            γ1  = gamma1(mat, σp[1], σp[2])
            fcm = γ1*mat.fc

            Ep1, Ep2, Ep3 = principal_isotropic_moduli(mat, εp, γ1)

            κ = 0.4

            if σp[3] >= κ*fcm # σc'  low compression
                println("Baixa compressão ou tração")
                println("D(isotrop_modif)")
                Et = ( abs(σp[1])*Ep1 + abs(σp[2])*Ep2 + abs(σp[3])*Ep3 ) / (abs(σp[1]) + abs(σp[2]) + abs(σp[3]))
                Ep1 = Ep2 = Ep3 = E12 = E23 = E13 = Et
            else
                println("Alta compressão")
                println("D(ortho_modif)")
                E12, E23, E13 = orthotropic_moduli(Ep1, Ep2, Ep3, σp)
            end

            # Orthotropic D matrix considering plane stress conditions
            Dp = 1.0/(1-ν^2)*[ Ep1    ν*E12   ν*E13  0.0        0.0        0.0 
                               ν*E12  Ep2     ν*E23  0.0        0.0        0.0 
                               ν*E13  ν*E23   Ep3    0.0        0.0        0.0 
                               0.0    0.0     0.0    (1-ν)*E12  0.0        0.0
                               0.0    0.0     0.0    0.0        (1-ν)*E13  0.0  
                               0.0    0.0     0.0    0.0        0.0        (1-ν)*E23 ]
        end

        # Fix D
        @show nfplanes
        for i in nfplanes
            if nfplanes == 0 # Dp não se modifica
                continue
            else
                @show ipd.fplanes
                active_fails = [ find_dir(V, p.V) for p in ipd.fplanes if p.active ]
                @show active_fails
                fixD(mat, Dp, active_fails)
            end
        end
        D = R'*Dp*R # rotates tensor Dp to xyz system
        return D
    end
end

function stress_update(mat::Orthotropic, ipd::OrthotropicIpState, Δε::Array{Float64,1})
    println("------------------------------------------")
    println("Stress Update")
    σ0 = copy(ipd.σ)
    De = calcDe(mat.E0, mat.ν, ipd.shared_data.model_type)
    σtr = ipd.σ + De*Δε

    E = mat.E0
    ν = mat.ν

    f = loading_func(mat, σtr)
    ipd.unloading = (f<ipd.fmax)
    ipd.fmax = max(f, ipd.fmax)

    nactive_fails = length(ipd.fplanes)

    nfplanes = length(ipd.fplanes)
    if nfplanes>0
        nactive_fails = sum( p.active for p in ipd.fplanes )
    else
        nactive_fails = 0
    end

    if !ipd.crushed && nactive_fails==0 # no fails
        println("Não está crushed e não tem planos de falha")
        σp, V = eigen(ipd.σ)
        R  = zeros(6,6)
        tensor_rot!(V, R)
        εp = R*ipd.ε

        if ipd.unloading
            println("Descarregando")
            println("De")
            D = calcDe(mat.E0, mat.ν, ipd.shared_data.model_type)
        else
            println("Carregando")
            p  = sortperm(σp, rev=true)
            σp = σp[p] # ordered stresses
            V  = V[:,p]

            R = zeros(6,6)
            tensor_rot!(V, R)
            εp  = R*ipd.ε # strain associated with σp

            γ1  = gamma1(mat, σp[1], σp[2])
            fcm = γ1*mat.fc

            # Use secant Young's moduli
            Δεp = R*Δε # incremental strain
            Ep1, Ep2, Ep3 = principal_orthotropic_moduli(mat, εp, Δεp, γ1)

            κ = 0.4
            if σp[3] >= κ*fcm # σc'  low compression
                println("Baixa compressão ou tração")
                println("D(isotrop)")
                if σp[1]!=0 || σp[2]!=0 || σp[3]!=0 
                    Et = ( abs(σp[1])*Ep1 + abs(σp[2])*Ep2 + abs(σp[3])*Ep3 ) / (abs(σp[1]) + abs(σp[2]) + abs(σp[3]))
                else
                    Et = mat.E0
                end
                D = calcDe(Et, mat.ν, ipd.shared_data.model_type)


            else # σp[3]<κ*fcmax high compression
                println("Alta compressão")
                println("D(ortho)")
                E12, E23, E13 = orthotropic_moduli(Ep1, Ep2, Ep3, σp)

                # Orthotropic D matrix
                # Notice that Amaru considers a general shear stress components (e.g. εxy)
                # and not the engineering definitions (e.g, γxy).
                Dp = 1/((1+ν)*(1-2*ν))*[ (1-ν)*Ep1   ν*E12       ν*E13      0.0          0.0          0.0 
                                         ν*E12       (1-ν)*Ep2   ν*E23      0.0          0.0          0.0 
                                         ν*E13       ν*E23      (1-ν)*Ep3   0.0          0.0          0.0 
                                         0.0         0.0         0.0        (1-2*ν)*E12  0.0          0.0
                                         0.0         0.0         0.0        0.0          (1-2*ν)*E13  0.0  
                                         0.0         0.0         0.0        0.0          0.0          (1-2*ν)*E23 ]
                D = R'*Dp*R # rotates tensor Dp to xyz system

            end

        end
        Δσ = D*Δε


    elseif ipd.crushed && Δε[3] < 0.0 ###Material failed already in COMPRESSION###
        println("Crushed")
        if nfplanes==0
            println("Planos de falha = 0")
            σp, V = eigen(ipd.σ)
        elseif nfplanes==1
            println("Plano de falha definido = 1")
            σp, V = eigen_with_fixed_dir(ipd.σ, ipd.fplanes[1].V) # V1 should be the first column of V
        else
            println("Planos de falha definido = 2 ou 3")
            V = [ ipd.fplanes[1].V ipd.fplanes[2].V ipd.fplanes[3].V ]
            R = zeros(6,6)
            tensor_rot!(V, R)
            σp = R*ipd.σ # stresses associated with V1, V2 and V3
            σp = σp[1:3]
        end

        p  = sortperm(σp, rev=true)
        σp = σp[p] # ordered stresses
        V  = V[:,p]

        R  = zeros(6,6)
        tensor_rot!(V, R)
        εp =  R*ipd.ε
        Δεp = R*Δε # incremental strain


        if minimum((εp+Δεp)[1:3]) < mat.εu # stress release
            println("Stress release")
            Δσ = -ipd.σ
        else
            if Δεp[3] > 0.0 # unloading (different than the original paper)
                println("unloading")
                D = calcDe(mat.E0, mat.ν, ipd.shared_data.model_type)
                println("Ee")
            else
                println("loading")
                if ipd.unloading && εp[3] > ipd.ep3max && σp[3] > ipd.sp3max
                #if ipd.unloading && εp[3] > ipd.ep3max
                #if ipd.unloading
                    D = calcDe(mat.E0, mat.ν, ipd.shared_data.model_type)
                    println("E(elastic)")
                else
                    Et = ( sigma(mat, εp[3]+Δεp[3]) - sigma(mat, εp[3]) ) / Δεp[3]
                    D  = calcDe(Et, mat.ν, ipd.shared_data.model_type)
                    println("Et")
                end
            end
            Δσ = D*Δε
        end

    else ### Material failed already in TENSION ###
        println("Tração")
        if nfplanes==0
            println("Planos de falha = 0")
            σp, V = eigen(ipd.σ)
        elseif nfplanes==1
            println("Plano de falha definido = 1")
            σp, V = eigen_with_fixed_dir(ipd.σ, ipd.fplanes[1].V ) # V1 should be the first column of V
        else
            println("Planos de falha definidos = 2 ou 3")
            V = [ ipd.fplanes[1].V ipd.fplanes[2].V ipd.fplanes[3].V ]
            R = zeros(6,6)
            tensor_rot!(V, R)
            σp = R*ipd.σ # stresses associated with V1, V2 and V3
            σp = σp[1:3]
        end


        # sort principal stresses
        p  = sortperm(σp, rev=true)
        σp = σp[p] # ordered stresses
        V  = V[:,p]
        R  = zeros(6,6)
        tensor_rot!(V, R)
        εp  = R*ipd.ε # strain associated with σp

        if ipd.unloading
            println("Descarregando")
            println("De")

            Dp = calcDe(mat.E0, mat.ν, ipd.shared_data.model_type)
        else # loading
            println("Carregando")

            γ1  = gamma1(mat, σp[1], σp[2])
            fcm = γ1*mat.fc

            # Use secant Young's moduli
            Δεp = R*Δε # incremental strain
            Ep1, Ep2, Ep3 = principal_orthotropic_moduli(mat, εp, Δεp, γ1)

            κ = 0.4

            if σp[3] >= κ*fcm # σc'  low compression
                println("Baixa compressão ou tração")
                println("D(isotrop_modif)")
                Et = ( abs(σp[1])*Ep1 + abs(σp[2])*Ep2 + abs(σp[3])*Ep3 ) / (abs(σp[1]) + abs(σp[2]) + abs(σp[3]))
                Ep1 = Ep2 = Ep3 = E12 = E23 = E13 = Et
            else  # σp[3]<κ*fcmax high compression
                println("Alta compressão")
                println("D(ortho_modif)")
                E12, E23, E13 = orthotropic_moduli(Ep1, Ep2, Ep3, σp)
            end

            # Orthotropic D matrix considering plane stress conditions
            Dp = 1.0/(1-ν^2)*[ Ep1    ν*E12   ν*E13  0.0        0.0        0.0 
                               ν*E12  Ep2     ν*E23  0.0        0.0        0.0 
                               ν*E13  ν*E23   Ep3    0.0        0.0        0.0 
                               0.0    0.0     0.0    (1-ν)*E12  0.0        0.0
                               0.0    0.0     0.0    0.0        (1-ν)*E13  0.0  
                               0.0    0.0     0.0    0.0        0.0        (1-ν)*E23 ]
        end

        # Fix D
        for i in nfplanes
            if nfplanes == 0 # Dp não se modifica
                continue
            else
                active_fails = [ find_dir(V, p.V) for p in ipd.fplanes ]
                fixD(mat, Dp, active_fails)
            end
        end
        D = R'*Dp*R # rotates tensor Dp to xyz system
        Δσ = D*Δε
    end

    ipd.ε += Δε
    ipd.σ += Δσ

    ep = minimum(εp[1:3])
    sp = minimum(σp[1:3])

    if ep < ipd.ep3max
        ipd.ep3max = ep
        ipd.sp3max = sp
    else 
        ipd.ep3max = ipd.ep3max
        ipd.sp3max = ipd.sp3max
    end

    
    # Check for new failure planes
    println("--------------------")
    println("Calculate a new σp and V")

    if nfplanes==0
        println("Planos de falha = '0'")
        σp, V = eigen(ipd.σ)
    elseif nfplanes==1
        println("Planos de falha = '1'")
        σp, V = eigen_with_fixed_dir(ipd.σ, ipd.fplanes[1].V ) # V1 should be the first column of V
    else # nfplanes == 3
        println("Planos de falha = '2 ou 3'")
        V = [ ipd.fplanes[1].V ipd.fplanes[2].V ipd.fplanes[3].V ]
    end
    R  = zeros(6,6)
    tensor_rot!(V, R)
    σp = R*ipd.σ # stresses associated with V1, V2 and V3
    εp = R*ipd.ε # strain associated with σp

    println("--------------------")
    println("Check for new failure planes")
    nfplanes = length(ipd.fplanes)
    for i=1:3
        if σp[i]>mat.ft # new fixed plane
            if nfplanes==0
                println("agora tem 1 plano fixo")
                plane = FixedPlane(V[:,i], true, true, εp[i])
                push!(ipd.fplanes, plane)
            elseif nfplanes==1
                if norm(V[:,i] - ipd.fplanes[1].V)<0.001 # V[:,i] == V1 plano já existe
                    continue
                else
                    println("agora tem 2 plano fixo")
                    plane = FixedPlane(V[:,i], true, true, εp[i])
                    push!(ipd.fplanes, plane)
                end
            elseif nfplanes==2
                if norm(V[:,i] - ipd.fplanes[1].V)<0.001 && norm(V[:,i] - ipd.fplanes[2].V)<0.001 
                    continue
                else
                    println("agora tem 3 plano fixo")
                    plane = FixedPlane(V[:,i], true, true, εp[i])
                    push!(ipd.fplanes, plane)
                end
            end
        end
    end

    # Check for active failure planes
    println("--------------------")
    println("Check for active failure planes")
    nfplanes = length(ipd.fplanes)
    if nfplanes>0
        for i=1:nfplanes
            plane = ipd.fplanes[i]
            dir = find_dir(V, plane.V)
            if εp[dir]>=plane.εf 
                plane.active = true
                σp[dir] = 0.0
            else
                plane.active = false
            end
        end

        ipd.σ = R'*σp
        Δσ = ipd.σ - σ0
    end
    
    # Check for compression crushing
    if ipd.crushed==false
        for i=1:3
            if σp[i]<mat.fc
                ipd.crushed = true
                break
            end
        end
    end


    return Δσ
end


function ip_state_vals(mat::Orthotropic, ipd::OrthotropicIpState)
    ndim  = ipd.shared_data.ndim
    σ, ε  = ipd.σ, ipd.ε
    j1    = tr(σ)
    srj2d = √J2D(σ)

    D = stress_strain_dict(σ, ε, ndim)

    return D
end
