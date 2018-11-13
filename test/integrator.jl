using Amaru

mat = Orthotropic(E=31.72e6, nu=0.2, fc=30.68e3, ft=3.1e3, epsc=0.00218, epsu= 0.00313, fu=26.48e3)

int = MechIntegrator(mat)

nu = 0.2
#Δε = [ 0.0025*nu, 0.0025*nu, -0.0025, 0, 0, 0 ]
#Δε = [ 0.0000, 0.0000, -0.0033, 0, 0, 0 ]*-0.025
#stress_update(int, Δε, nincs=20)

Δε = [ 0.0000, 0.0000, -0.0033, 0, 0, 0 ]*0.85
stress_update(int, Δε, nincs=30)

Δε = [ 0.0000, 0.0000, -0.0033, 0, 0, 0 ]*-0.450
stress_update(int, Δε, nincs=30)

Δε = [ 0.0000, 0.0000, -0.0033, 0, 0, 0 ]*0.55
stress_update(int, Δε, nincs=30)

using PyPlot
t = int.table
grid("on", linewidth=0.5, color="lightgray", linestyle="-")     # propriedades do grid
plot(t[:ezz], t[:szz], "-o")
tight_layout()                                                  # centrar a figura
#@show t
#;
