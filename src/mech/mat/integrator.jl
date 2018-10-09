export MechIntegrator
export stress_update

mutable struct MechIntegrator
    ipd::IpState
    mat::Material
    table::DTable
    function MechIntegrator(mat::Material)
        this = new()
        shared_data = SharedAnalysisData()
        this.ipd = new_ip_state(mat, shared_data)
        this.mat = mat
        this.table = DTable()
        return this
    end
end


function stress_update(int::MechIntegrator, Δε; nincs=1)
    Δεi = Δε/nincs
    for i=1:nincs
        stress_update(int.mat, int.ipd, Δεi)
        vals = ip_state_vals(int.mat, int.ipd)
        push!(int.table, vals)
    end
end

