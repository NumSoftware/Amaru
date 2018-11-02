# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct MechSolid<:Mechanical
    id    ::Int
    shape ::ShapeType
    cell  ::Cell
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::TagType
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    shared_data::SharedAnalysisData

    function MechSolid(); 
        return new() 
    end
end

matching_shape_family(::Type{MechSolid}) = SOLID_SHAPE

function elem_init(elem::MechSolid)
    ipdata_ty = typeof(elem.ips[1].data)
    if :h in fieldnames(ipdata_ty)
        # Element volume/area
        V = 0.0
        C = elem_coords(elem)
        for ip in elem.ips
            dNdR = elem.shape.deriv(ip.R)
            J    = dNdR*C
            detJ = det(J)
            V   += detJ*ip.w
        end

        # Representative length size for an integration point
        nips = length(elem.ips)
        ndim = elem.shared_data.ndim
        h = (V/nips)^(1/ndim)

        for ip in elem.ips
            ip.data.h = h
        end
    end

    return nothing
end

function distributed_bc(elem::MechSolid, facet::Union{Facet, Nothing}, key::Symbol, fun::Functor)
    ndim  = elem.shared_data.ndim

    # Check bcs
    (key == :tz && ndim==2) && error("distributed_bc: boundary condition $key is not applicable in a 2D analysis")
    !(key in (:tx, :ty, :tz, :tn)) && error("distributed_bc: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem))")

    target = facet!=nothing ? facet : elem
    nodes  = target.nodes
    nnodes = length(nodes)
    t      = elem.shared_data.t

    # Force boundary condition
    nnodes = length(nodes)

    # Calculate the target coordinates matrix
    C = nodes_coords(nodes, ndim)

    # Vector with values to apply
    Q = zeros(ndim)

    # Calculate the nodal values
    F     = zeros(nnodes, ndim)
    shape = target.shape
    ips   = get_ip_coords(shape)

    for i=1:size(ips,1)
        R = vec(ips[i,:])
        w = R[end]
        N = shape.func(R)
        D = shape.deriv(R)
        J = D*C
        nJ = norm2(J)
        X = C'*N
        if ndim==2
            x, y = X
            val = fun(t,x,y,0.0)
            if key == :tx
                Q = [val, 0.0]
            elseif key == :ty
                Q = [0.0, val]
            elseif key == :tn
                n = [J[1,2], -J[1,1]]
                Q = val*n/norm(n)
            end
        else
            x, y, z = X
            val = fun(t,x,y,z)
            if key == :tx
                Q = [val, 0.0, 0.0]
            elseif key == :ty
                Q = [0.0, val, 0.0]
            elseif key == :tz
                Q = [0.0, 0.0, val]
            elseif key == :tn && ndim==3
                n = cross(J[1,:], J[2,:])
                Q = val*n/norm(n)
            end
        end
        F += N*Q'*(nJ*w) # F is a matrix
    end

    # generate a map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in target.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end

function setB(shared_data::SharedAnalysisData, dNdX::Matx, detJ::Float64, B::Matx)
    ndim, nnodes = size(dNdX)
    B   .= 0.0

    if ndim==2
        for i in 1:nnodes
            j = i-1
            B[1,1+j*ndim] = dNdX[1,i]
            B[2,2+j*ndim] = dNdX[2,i]
            B[6,1+j*ndim] = dNdX[2,i]/SR2; B[6,2+j*ndim] = dNdX[1,i]/SR2
        end
        if shared_data.model_type==:axisymmetric
            for i in 1:nnodes
                N =elem.shape.func(R)
                j = i-1
                r = R[0]
                B[1,1+j*ndim] = dNdX[1,i]
                B[2,2+j*ndim] = dNdX[2,i]
                B[3,1+j*ndim] =    N[i]/r
                B[6,1+j*ndim] = dNdX[2,i]/SR2; B[6,2+j*ndim] = dNdX[1,i]/SR2
            end
        end
    else
        for i in 1:nnodes
            dNdx = dNdX[1,i]
            dNdy = dNdX[2,i]
            dNdz = dNdX[3,i]
            j    = i-1
            B[1,1+j*ndim] = dNdx
            B[2,2+j*ndim] = dNdy
            B[3,3+j*ndim] = dNdz
            B[4,2+j*ndim] = dNdz/SR2;   B[4,3+j*ndim] = dNdy/SR2
            B[5,1+j*ndim] = dNdz/SR2;   B[5,3+j*ndim] = dNdx/SR2
            B[6,1+j*ndim] = dNdy/SR2;   B[6,2+j*ndim] = dNdx/SR2
        end
    end

    return detJ
end

function elem_stiffness(elem::MechSolid)
    ndim   = elem.shared_data.ndim
    th     = elem.shared_data.thickness
    nnodes = length(elem.nodes)
    C = elem_coords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(6, nnodes*ndim)

    DB = Array{Float64}(undef, 6, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)

    for ip in elem.ips

        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        @gemm dNdX = inv(J)*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        setB(elem.shared_data, dNdX, detJ, B)

        # compute K
        coef = detJ*ip.w*th
        D    = calcD(elem.mat, ip.data) 
        @gemm DB = D*B
        @gemm K += coef*B'*DB
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end

function elem_mass(elem::MechSolid)
    ndim   = elem.shared_data.ndim
    th     = elem.shared_data.thickness
    nnodes = length(elem.nodes)
    ρ = elem.mat.ρ
    C = elem_coords(elem)
    M = zeros(nnodes*ndim, nnodes*ndim)
    N = zeros(ndim, nnodes*ndim)
    J = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        # compute N matrix
        Ni   = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)

        for i=1:nnodes
            for j=1:ndim
                N[j, (i-1)*ndim+j] = Ni[i]
            end
        end

        @gemm J = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute K
        coef = ρ*detJ*ip.w*th
        @gemm M += coef*N'*N
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    return M, map, map
end


function elem_update!(elem::MechSolid, U::Array{Float64,1}, F::Array{Float64,1}, Δt::Float64)
    ndim   = elem.shared_data.ndim
    th     = elem.shared_data.thickness
    nnodes = length(elem.nodes)
    keys   = (:ux, :uy, :uz)[1:ndim]
    map    = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    dU = U[map]
    dF = zeros(nnodes*ndim)
    B  = zeros(6, nnodes*ndim)

    DB = Array{Float64}(undef, 6, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)
    Δε = zeros(6)

    C = elem_coords(elem)
    for ip in elem.ips

        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        @gemm dNdX = inv(J)*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        setB(elem.shared_data, dNdX, detJ, B)

        @gemv Δε = B*dU
        Δσ   = stress_update(elem.mat, ip.data, Tensor2(Δε))
        coef = detJ*ip.w*th
        dF += coef*B'*Δσ
    end

    F[map] += dF
end

