using Amaru

mat = Orthotropic(E=31.72e6, nu=0.2, fc=30.68e3, ft=3.1e3, epsc=0.00218, epsu= 0.00313, fu=26.48e3)

int = MechIntegrator(mat)

nu = 0.2
euu= -0.0033

Δε = [ euu*nu, euu*nu, euu, 0, 0, 0 ]*0.85
stress_update(int, Δε, nincs=100)

Δε = [ euu*nu, euu*nu, euu, 0, 0, 0 ]*-0.45
stress_update(int, Δε, nincs=100)

Δε = [ euu*nu, euu*nu, euu, 0, 0, 0 ]*0.55
stress_update(int, Δε, nincs=100)

using PyPlot
t = int.table
grid("on", linewidth=0.5, color="lightgray", linestyle="-")     # propriedades do grid
plot(t[:ezz], t[:szz], "-o")
tight_layout()                                                  # centrar a figura

