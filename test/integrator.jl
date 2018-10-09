using Amaru

mat = Orthotropic(E=31.72e6, nu=0.2, fc=30.68e3, ft=3.1e3, epsc=0.00218, epsu= 0.00313, fu=26.48e3)

int = MechIntegrator(mat)

nu = 0.2
Δε = [ 0.0025*nu, 0.0025*nu, -0.0025, 0, 0, 0 ]

stress_update(int, Δε, nincs=40)

using PyPlot
t = int.table
plot(t[:ezz], t[:szz], "-o")

@show t
