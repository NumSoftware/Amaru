# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

include("elem/thermomech.jl")
export elem_config_dofs, elem_init, elem_stiffness, elem_update!, elem_vals
export set_state

# Thermomechanical Elements
include("elem/thermomech-solid.jl")
include("elem/thermo-solid.jl")
include("elem/thermomech-shell.jl")

# Models for solid elements (2D and 3D)
#include("mat/elastic-solid-lin-cond.jl")
include("mat/lin-thermo.jl")
include("mat/elastic-solid-thermo.jl")
include("mat/elastic-shell-thermo.jl")

include("thermomech-solver.jl")

export tm_solve!
