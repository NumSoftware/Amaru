# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct SharedAnalysisData
    ndim::Int              # Analysis dimension
    model_type::Symbol     # Analysis type (eg. :general, :plane_stress, :plane_strain, etc.
    thickness::Float64     # Model thickness
    t::Float64             # Time in time dependent analysis
    function SharedAnalysisData()
        this = new()
        this.ndim = 3
        this.model_type = :general
        this.t = 0.0
        this.thickness = 1.0
        return this
    end
end

