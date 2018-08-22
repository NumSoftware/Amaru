# From the book: 
# Finite Elements for Structural Analysis. pg 168
# William Weaver & Paul Johnston

using Amaru
using Base.Test

# Mesh generation

block = Block3D( [0 0 0; 1 1 1], nx=1, ny=1, nz=1, shape=HEX8) 
mesh = Mesh(block, verbose=true, reorder=false)

# Domain definition

materials = [
    MaterialBind(0, ElasticSolid(E=1.0, nu=0.3) ),
]

# Load cases

bcs1 = [
    NodeBC(:(x==0 && y==0 && z==0), :(ux=0, uy=0) ),
    NodeBC(:(x==1 && y==0 && z==0), :(uy=0) ),
    NodeBC(:(x==0 && y==1 && z==0), :(ux=0) ),
    NodeBC(:(z==0), :(uz=0) ),
    NodeBC(:(z==1), :(fz=1) ),
]

bcs2 = [
    NodeBC(:(x==0 && y==0 && z==0), :(ux=0, uy=0) ),
    NodeBC(:(x==1 && y==0 && z==0), :(uy=0) ),
    NodeBC(:(x==0 && y==1 && z==0), :(ux=0) ),
    NodeBC(:(z==0), :(uz=0) ),
    EdgeBC(:(y==1 && z==1), :(ty=2) ),
]

bcs3 = [
    NodeBC(:(x==0 && y==0 && z==0), :(ux=0, uy=0) ),
    NodeBC(:(x==1 && y==0 && z==0), :(uy=0) ),
    NodeBC(:(x==0 && y==1 && z==0), :(ux=0) ),
    NodeBC(:(z==0), :(uz=0) ),
    FaceBC(:(x==1), :(tx=3*z) ),
]

bcs4 = [
    NodeBC(:(x==0 && y==0 && z==0), :(ux=0, uy=0) ),
    NodeBC(:(x==1 && y==0 && z==0), :(uy=0) ),
    NodeBC(:(x==0 && y==1 && z==0), :(ux=0) ),
    NodeBC(:(z==0), :(uz=0) ),
    ElemBC(:(x>=0), :(tz=-1) ),
]

ana_list = ["Nodal load", "Edge load", "Triangular face load", "Volume load"]
bcs_list = [bcs1, bcs2, bcs3, bcs4]
dis_list = [
            [ 0.0, 0.0, 0.0, 0.0, 4.0    ,    4.0,     4.0,      4.0 ],
            [ 0.0, 0.0, 0.0, 0.0, 3.32088, 3.1998, -4.6002, -4.32047 ],
            [ 0.0, 0.0, 0.0, 0.0, 1.51044,-2.4501,  1.4499, -2.31023 ],
            [ 0.0, 0.0, 0.0, 0.0, -0.5   ,   -0.5,    -0.5,     -0.5 ] ]

for (ana, bcs, dis) in zip(ana_list, bcs_list, dis_list)

    println("\nLoad case: $ana \n")

    dom = Domain(mesh, materials)
    solve!(dom, bcs, nincs=1, nouts=1, verbose=true)

    println("Displacements:")
    D = nodes_dof_vals(dom.nodes)[[:ux, :uy, :uz]]
    println(D)
    println(D[:uz])
    println(dis)

    @test dis ≈ D[:uz] atol=1e-5

    println("Stress:")
    S = elems_ip_vals(dom.elems[1])[[:sxx, :syy, :szz, :syz, :sxz, :sxy]]
    println(S)

    println("Support reactions:")
    F = nodes_dof_vals(dom.nodes[:(z==0)])[[:fx, :fy, :fz]]
    println(F)
end


