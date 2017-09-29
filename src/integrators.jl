struct NLSEIntegrator end

# placeholder for future cache-ing
function initialize!(integrator)

end

@muladd function perform_step!(integrator, repeat_step=false)
    @unpack t, dt, z, dz, uprev, u, f = integrator

end
