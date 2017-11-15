function get_nonlinearcoeff(res, mode, source)
    getω(source)/c * res.material.nonlinearindex * mode.corefraction(source) / mode.effectivearea(source)
end

function get_ωmesh(tmesh)
    dt_mesh = tmesh[2]-tmesh[1]
    num_pts = length(tmesh)
    ω = collect(2pi*(-num_pts/2:num_pts/2-1)/(num_pts*dt_mesh))
    ω
end
