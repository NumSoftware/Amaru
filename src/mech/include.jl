# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

include("elem/mech.jl")
export elem_config_dofs, elem_init, elem_stiffness, elem_mass, elem_update!, elem_vals
export set_state

# Elements
include("elem/mech-solid.jl")
#include("elem/in-mech-solid.jl")
include("elem/mech-rod.jl")
include("elem/mech-embrod.jl")
include("elem/mech-beam.jl")
include("elem/mech-joint.jl")
include("elem/mech-joint1d.jl")


# Models for solid elements (2D and 3D)
include("mat/elastic-solid.jl")
#include("mat/in-elastic-solid.jl")
include("mat/dp-solid.jl")
include("mat/vm-solid.jl")
include("mat/mazars-solid.jl")
include("mat/sc-solid.jl")
include("mat/ortho-solid.jl")

# Models for truss elements
include("mat/elastic-rod.jl")
include("mat/pp-rod.jl")

# Models for beams
include("mat/elastic-beam.jl")

# Models for joint elements
include("mat/elastic-joint.jl")
include("mat/mc-joint.jl")

# Models for 1D joint elements
include("mat/elastic-joint1d.jl")
include("mat/ceb-joint1d.jl")

# Stress-update integrator
include("mat/integrator.jl")

# Solvers
include("solver.jl")
include("dyn-solver.jl")
export solve!
