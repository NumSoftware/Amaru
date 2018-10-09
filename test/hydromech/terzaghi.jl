using Amaru
using Test


# Mesh generation
# ===============

blocks = [
    Block3D( [0 0 0; 1 1 10], nx=1, ny=1, nz=10, cellshape=HEX8),
]

msh = Mesh(blocks, verbose=true)


# Finite element analysis
# =======================

# Analysis data
load = 10.0   
k    = 1.0E-5  # permeability
E    = 5000.0  # Young modulus
nu   = 0.25    # Poisson
gw   = 10.0    # water specific weight
hd   = 10.0    # drainage height
mv   = (1+nu)*(1-2*nu)/(E*(1-nu))
cv   = k/(mv*gw)   # consolidation coefficient
T = [ 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0 ]
times = T*(hd^2/cv)
t1 = times[1]/10
times .+= t1

materials = [
    MaterialBind(:solids, ElasticSolidLinSeep(E=E, nu=nu, k=k, gw=gw) ),
]

logger = [
    NodeGroupLogger(:(x==0 && y==0) ),
]

dom = Domain(msh, materials, logger)


# Stage 1: hydrostatic pore-pressure
# ==================================

tlong = 10000*hd^2/cv

bcs = [
    NodeBC(:(z==0), :(ux=0, uy=0, uz=0) ),
    NodeBC(:(x>=0), :(ux=0, uy=0) ),
    NodeBC(:(z==10), :(uw=0.) ),
]

hm_solve!(dom, bcs, end_time=tlong, nincs=2, tol=1e-2, nouts=1, verbose=true)

dom.shared_data.t = 0.0


# Stage 2: loading
# ================

pt(t) = -load
#pt(t) = t>=t1? -load : -load/t1*t

bcs = [
    NodeBC(:(z==0), :(ux=0, uy=0, uz=0) ),
    NodeBC(:(x>=0), :(ux=0, uy=0) ),
    FaceBC(:(z==10), :(tz=$pt(t)) ),
    NodeBC(:(z==10), :(uw=0.) ),
]

hm_solve!(dom, bcs, end_time=t1, nincs=4, tol=1e-2, nouts=1, verbose=true)


# Stage 3: draining
# =================

bcs = [
    NodeBC(:(z==0), :(ux=0, uy=0, uz=0) ),
    NodeBC(:(x>=0), :(ux=0, uy=0) ),
    FaceBC(:(z==10), :(tz=$pt(t)) ),
    NodeBC(:(z==10), :(uw=0.) ),
]

Uw_vals = [] # A list with porepressure vectors
for t in times[1:end]
    hm_solve!(dom, bcs, end_time=t, nincs=20, tol=1e-2, nouts=1, verbose=false)
    push!( Uw_vals, logger[1].book[end][:uw] )
end


# Output
# ======

if Amaru.Debug.makeplots
    # Terzaghi's 1-d consolidation
    function calc_Ue(Z, T)
        sum = 0.0
        for i=0:4
            M = pi/2*(2*i+1)
            sum = sum + 2/M*sin(M*Z)*exp(-M^2*T)
        end
        return sum
    end

    using PyPlot


    # numerical curves
    book  = logger[1].book
    Uwini = book.tables[2][:uw]       # hydrostatic porepressure
    Z     = 1 .- book.tables[2][:z]/hd # normalized depth
    for Uw in Uw_vals
        dUw = (Uw - Uwini)/load       # excess of porepressure
        plot(dUw, Z, "-o")
    end

    # analytical curves
    for Ti in T
        Ue = calc_Ue.(Z, Ti)
        plot(Ue, Z, "k")
    end

    ylim(1,0)
    show()
end
