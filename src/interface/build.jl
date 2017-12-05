function build_problem(model, probtype::DynamicNLSE;
                       time_pts=2^13, time_window=10.0ps, kwargs...)

    tmesh = linspace(-time_window/2, time_window/2, time_pts)
    zspan = (0.0, pathlength(model))

    ω_zeroed = get_ωmesh(tmesh)
    const dt_mesh = tmesh[2]-tmesh[1]
    ω_zeroed = fftshift(ω_zeroed)

    const ω0 = get_ω0(model)
    const ω = ω_zeroed + ω0
    u0 = derive_pulse(model, tmesh)

    f, g, D, planned_fft!, planned_ifft! = build_GNLSE(model, ω, tmesh)
    planned_ifft! * view(u0, :, 1)
    prob = ODEProblem(
        f,
        u0,
        zspan
    )
    prob_exp = DynamicNLSEProblem(
        prob,
        model,
        fftshift(ω),
        ω0,
        tmesh,
        planned_fft!,
        planned_ifft!,
        D
    )
    return prob_exp
end

function build_problem(model, probtype::StochasticNLSE;
                       time_pts=2^13, time_window=10.0ps, kwargs...)

    tmesh = linspace(-time_window/2, time_window/2, time_pts)
    zspan = (0.0m, pathlength(model))

    ω_zeroed = get_ωmesh(tmesh)
    const dt_mesh = tmesh[2]-tmesh[1]
    ω_zeroed = fftshift(ω_zeroed)

    const ω0 = frequency(model)
    const ω = ω_zeroed + ω0
    u0 = derive_pulse(model, tmesh)
    pulse_unit = oneunit(first(u0))
    length_unit = oneunit(first(zspan))
    u0 /= pulse_unit

    f, g, D, planned_fft!, planned_ifft! = build_GNLSE(model, ω, tmesh)
    planned_ifft! * view(u0, :, 1)
    prob = SDEProblem(
        f,
        g,
        u0,
        (zspan[1]/length_unit, zspan[2]/length_unit)
    )
    prob_exp = StochasticNLSEProblem(
        prob,
        model,
        fftshift(ω),
        ω0,
        tmesh,
        planned_fft!,
        planned_ifft!,
        D
    )
    return prob_exp
end

function build_problem(model, probtype::DynamicIkeda;
                       time_pts=2^13, time_window=10.0, kwargs...)

    tmesh = linspace(-time_window/2, time_window/2, time_pts)
    zspan = (0.0m, pathlength(model))

    ω_zeroed = get_ωmesh(tmesh)
    const dt_mesh = tmesh[2]-tmesh[1]
    ω_zeroed = fftshift(ω_zeroed)
    const ω0 = model.ω0
    const ω = ω0 + ω_zeroed

    u0 = derive_pulse(model, tmesh)

    f, D, planned_fft!, planned_ifft! = build_GNLSE(model, ω, tmesh)
    ikeda_callback = build_ikedacallback(model)
    prob = ODEProblem(f,
                      planned_ifft! * u0,
                      zspan)
    prob_exp = DynamicIkedaProblem(prob,
                                   model,
                                   fftshift(ω),
                                   ω0,
                                   tmesh,
                                   planned_fft!,
                                   planned_ifft!,
                                   D,
                                   ikeda_callback)
    return prob_exp
end


# takes a scalar and turns it into a callable
function get_func(attr::Number)
    f(x...) = attr
    return f
end
get_func(attr) = attr

function build_problem(model::ToyModel, ::DynamicLL;
                       time_pts=2^10, time_window=10.0, kwargs...)
    FSR = model.FSR
    tmesh = linspace(-1/2/FSR, 1/2/FSR, time_pts)
    τspan = (0.0, time_window)

    ω = get_ωmesh(tmesh)
    ω = fftshift(ω)
    ω0 = 0.0#toy model ω0 is 0
    ω .+= ω0

    #TODO: random start conditions don't seem to do anything
    u0 = derive_pulse(model.power_in, model.pulsetime).(tmesh)

    f, D, planned_fft!, planned_ifft! = build_GLLE(model, ω, tmesh)

    prob = ODEProblem(f, planned_ifft! * u0, τspan; kwargs...)
    prob_LL = DynamicLLProblem(prob, model, fftshift(ω), ω0, tmesh, planned_fft!, planned_ifft!, D)
end

function build_GNLSE(model::ToyModel, ω, tmesh)

    const α = linearloss(model, ω)
    const γnl = nonlinearcoeff(model, ω)
    const beta = stationarydispersion(model, ω)

    const planned_fft! = plan_fft!(Vector((1+0im)*tmesh), flags=FFTW.MEASURE)
    const planned_ifft! = plan_ifft!(Vector((1+0im)*tmesh), flags=FFTW.MEASURE)
    const D = 1im .* beta .- α/2

    const use_raman = has_raman(model)

    const shock_term = build_model_shock(model, ω)
    const (raman_response, raman_noise) = build_model_raman(model, tmesh, ω)

    const frac_raman = 0.18

    function f(z, u, du)
        @. du = u * exp(D * z)
        planned_fft! * du
        if use_raman
            raman_term = planned_ifft! * abs2.(du) # prepare intensity for convolution
            # mulitiply in fourier space for convolution
            @. raman_term = raman_term * raman_response
            planned_fft! * raman_term # fourier transform back
            @. du = du * ((1-frac_raman)*abs2(du) + frac_raman*raman_term)
        else
            @. du = du * abs2(du)
        end
        planned_ifft! * du
        @. du = 1im * γnl * shock_term * du * exp(-D * z)
    end
    function g(z, u, du2) end
    return f, g, D, planned_fft!, planned_ifft!
end

function build_GNLSE(model::Model, ω, tmesh)
    # variables indexed to i.e. α[i,j,k] -> α[(ω), (modes), (structure segment)]
    const α = linearloss(model, ω)
    const γnl = nonlinearcoeff(model, ω)
    const beta = stationarydispersion(model, ω)

    unit_time = oneunit(1.0ps)
    fft_mesh = Vector(Complex.(tmesh./unit_time))
    const planned_fft! = plan_fft!(fft_mesh, flags=FFTW.MEASURE)
    const planned_ifft! = plan_ifft!(fft_mesh, flags=FFTW.MEASURE)
    const D = 1im .* beta .- α/2

    const use_raman = has_raman(model)

    const shock_term = build_model_shock(model, ω)
    const (raman_response, raman_noise) = build_model_raman(model, tmesh, ω)

    const num_mode = num_modes(model)
    function f(z, u, du)
        s_idx, s_curr   = currentstructure(model, z*m)
        raman_tens      = raman_tensor(s_curr, z*m)
        electronic_tens = electronic_tensor(s_curr, z*m)
        frac_raman      = raman_fraction(s_curr)
        med_curr        = s_curr.medium
        # First we change out of the interaction picture, turning all of our modes
        # from an interaction picture fourier-domain represenation to a time-domain
        # representation
        ut = zero(du)
        for mode_idx in 1:num_mode
            ut[:, mode_idx] .= view(u, :, mode_idx) .* exp.(view(D, :, mode_idx, s_idx) * z * m)
            planned_fft! * view(ut, :, mode_idx)
        end
        # Now we handle interactions between the modes in the time domain. Every
        # time du shows up, we are looking at the time-domain representaiton prepared
        # previously
        for mode_idx in 1:num_mode
            if use_raman[mode_idx]
                # We construct the raman noise, which is diagonal in the
                # frequency domain but multiplicative in our diff eq in the time domain
                if interacts_with_other_mode(med_curr.modes[mode_idx])
                    mode_pair_idx           = get_paired_mode_idx(med_curr, mode_idx)
                    mode_pol                = get_polarization(med_curr.modes[mode_idx])
                    mode_pair_pol           = get_polarization(med_curr.modes[mode_pair_idx])
                    overlap                 = med_curr.overlap[(mode_idx, mode_pair_idx)]
                    # mulitiply in fourier space for convolution
                    raman_straight_contrib  = raman_tens[mode_pol, mode_pol, mode_pol, mode_pol] * overlap * frac_raman
                    kerr_straight_contrib   = electronic_tens[mode_pol, mode_pol, mode_pol, mode_pol] * overlap * (1-frac_raman)
                    raman_cross_contrib     = raman_tens[mode_pol, mode_pair_pol, mode_pair_pol, mode_pol] * overlap * frac_raman
                    kerr_cross_contrib      = electronic_tens[mode_pol, mode_pair_pol, mode_pair_pol, mode_pol] * overlap * (1-frac_raman)
                    raman_straight_term     = planned_ifft! * abs2.(ut[:, mode_idx]) # prepare intensity for convolution - note that this allocates, so we don't call copy()
                    raman_cross_term        = planned_ifft! * (conj.(ut[:, mode_pair_idx]).*(ut[:, mode_idx])) # prepare intensity for convolution
                    raman_straight_term    .= raman_straight_term .* view(raman_response, :, mode_idx, s_idx)
                    planned_fft! * raman_straight_term # fourier transform back
                    raman_cross_term       .= raman_cross_term .* view(raman_response, :, mode_idx, s_idx)
                    planned_fft! * raman_cross_term # fourier transform back
                    # Raman noise is multiplicative in the time domain
                    du[:, mode_idx]        .= (view(ut, :, mode_idx) .* (kerr_straight_contrib*abs2.(view(ut, :, mode_idx)) + raman_straight_contrib.*raman_straight_term)
                                             + view(ut, :, mode_pair_idx) .* (kerr_cross_contrib.*conj.(view(ut, :, mode_pair_idx)).*(view(ut, :, mode_idx)) .+ raman_cross_contrib.*raman_cross_term))
                else
                    raman_term              = planned_ifft! * abs2.(ut[:, mode_idx]) # prepare intensity for convolution
                    # mulitiply in fourier space for convolution
                    raman_term             .= raman_term .* view(raman_response, :, mode_idx, s_idx)
                    planned_fft! * raman_term # fourier transform back
                    du[:, mode_idx]        .= view(ut, :, mode_idx) .* ((1-frac_raman)*abs2.(view(ut, :, mode_idx)) .+ frac_raman*(raman_term))
                end
            else
                du[:, mode_idx] .= view(du, :, mode_idx) .* abs2.(view(du, :, mode_idx))
            end
            planned_ifft! * view(du, :, mode_idx)
            du[:, mode_idx] .= (view(shock_term, :, mode_idx, s_idx)
                                  .* (1im * view(γnl, :, mode_idx, s_idx) .* view(du, :, mode_idx) * W * m)
                                 .* exp.(-view(D, :, mode_idx, s_idx) * z * m)
                                )
        end
    end

    function g(z, u, du2)
        s_idx, s_curr   = currentstructure(model, z*m)
        raman_tens      = raman_tensor(s_curr, z*m)
        electronic_tens = electronic_tensor(s_curr, z*m)
        frac_raman      = raman_fraction(s_curr)
        noise_unit      = oneunit(first(raman_noise))
        med_curr        = s_curr.medium
        # First we change out of the interaction picture, turning all of our modes
        # from an interaction picture fourier-domain represenation to a time-domain
        # representation
        ut = zero(du2)
        for mode_idx in 1:num_mode
            ut[:, mode_idx] .= view(u, :, mode_idx) .* exp.(view(D, :, mode_idx, s_idx) * z * m)
            planned_fft! * view(ut, :, mode_idx)
        end
        # Now we handle interactions between the modes in the time domain. Every
        # time du shows up, we are looking at the time-domain representaiton prepared
        # previously
        for mode_idx in 1:num_mode
            if use_raman[mode_idx]
                # TODO: in-place raman noise
                # We construct the raman noise, which is diagonal in the
                # frequency domain but multiplicative in our diff eq in the time domain
                raman_noise_term = fftshift(planned_fft! * (raman_noise[:, mode_idx, s_idx]/noise_unit))
                if interacts_with_other_mode(med_curr.modes[mode_idx])
                    mode_pair_idx           = get_paired_mode_idx(med_curr, mode_idx)
                    mode_pol                = get_polarization(med_curr.modes[mode_idx])
                    mode_pair_pol           = get_polarization(med_curr.modes[mode_pair_idx])
                    overlap                 = med_curr.overlap[(mode_idx, mode_pair_idx)]
                    # mulitiply in time domain
                    raman_straight_contrib  = raman_tens[mode_pol, mode_pol, mode_pol, mode_pol] * overlap * frac_raman
                    raman_cross_contrib     = raman_tens[mode_pol, mode_pair_pol, mode_pair_pol, mode_pol] * overlap * frac_raman
                    # Raman noise is multiplicative in the time domain
                    raman_noise_term      .*= view(ut, :, mode_idx) * raman_straight_contrib .+ view(ut, :, mode_pair_idx) * raman_cross_contrib
                else
                    # mulitiply in time domain
                    raman_noise_term      .*= view(ut, :, mode_idx) * frac_raman
                end
            end
            du2[:, mode_idx] = -(view(shock_term, :, mode_idx, s_idx)
                                 .* (planned_ifft! * raman_noise_term)
                                 .* exp.(-view(D, :, mode_idx, s_idx) * z * m)
                                )
        end
    end
    return f, g, D, planned_fft!, planned_ifft!
end

function build_GLLE(model, ω, tmesh)
    FSR = model.FSR
    α = linearloss(model, ω)
    γnl = nonlinearcoeff(model, ω)
    L = model.length
    detuning = model.detuning
    sqrtcoupling = sqrt(model.coupling)
    Ein = sqrt(model.power_in)
    beta = model.betacoeff
    beta_coeff = beta .* [(1im)^n/factorial(n) for n=0:(length(beta)-1)]
    const has_raman = model.has_raman
    const has_shock = model.has_shock

    const shock_term = build_model_shock(model, ω)
    const frac_raman = 0.18
    const (raman_response, raman_noise) = build_model_raman(model, tmesh, ω)

    const planned_fft! = plan_fft!(Vector((1+0im)*tmesh), flags=FFTW.MEASURE)
    const planned_ifft! = plan_ifft!(Vector((1+0im)*tmesh), flags=FFTW.MEASURE)
    const D = FSR * (-1im * L * Poly(beta_coeff, :ω).(ω) - α - 1im * detuning)
    #TODO: Macro for adding terms to function
    function f(z, u, du)
        @. du = (u * exp(D * z))
        planned_fft! * du
        if has_raman
            raman_term = planned_ifft! * copy(abs2.(du)) # prepare intensity for convolution
            # mulitiply in fourier space for convolution
            @. raman_term = raman_term * raman_response
            planned_fft! * raman_term # fourier transform back
            @. du = du * ((1-frac_raman)*abs2(du) + frac_raman*raman_term)
        else
            @. du = du * abs2(du)
        end
        planned_ifft! * du
        @. du = (1im * γnl * L * FSR * du) * exp(-D * z)
        du[1] = sqrtcoupling * Ein * FSR * exp(-D[1] * z)
    end
    return f, D, planned_fft!, planned_ifft!
end

function build_ikedacallback(model)
    L = model.length
    FSR = model.FSR
    const Ein = sqrt(model.power_in)
    const coupling = sqrt(model.coupling)
    const transmission = sqrt(1-coupling^2)
    const phase_shift = 2pi*L-model.detuning
    # summation is in fourier space
    function affect!(integrator)
        integrator.u = transmission * exp(1im*phase_shift) * integrator.u
        integrator.u[1] += coupling * Ein
    end
    ikeda_callback = PeriodicCallback(affect!, L)
end

function build_model_shock(model::ToyModel, ω)
    ω0 = model.ω0
    if model.has_shock
        shock_term = one(ω0) + ω/ω0
    else
        shock_term = one(ω0)
    end
    return shock_term
end
function build_model_shock(model::Model, ω)
    const ω0 = frequency(model.laser)
    shock_term = fill(one(ω0), length(ω), num_modes(model), num_structures(model))
    for (i, structure) in enumerate(model.waveguide.structures)
        for (j, mode) in enumerate(structure.medium.modes)
            if mode.has_shock
                shock_term[:, j, i] .+= ω/ω0
            end
        end
    end
    return shock_term
end

"""
    build_model_raman(model, tmesh, ω) -> (raman_response, raman_noise)

    Given a model and time points (tmesh), returns the calculated raman
    response for the material and its raman noise. Raman noise is diagonal
    in frequency domain.
"""
function build_model_raman(model::ToyModel, tmesh, ω)
    raman_response = (1+0im)*similar(tmesh)
    if model.has_raman
        tau1 = 0.0122; tau2 = 0.032
        raman_timeresponse = @. (1+0im)*(tmesh > 0) * (tau1^2 + tau2^2)/tau1/tau2^2*exp(-tmesh/tau2)*sin(tmesh/tau1)*length(tmesh)
        raman_response = ifft!(fftshift(raman_timeresponse))
    end
    return raman_response, zero(raman_response)
end
function build_model_raman(model::Model, tmesh, ω)
    ω0 = frequency(model)
    raman_response = zeros(Complex, length(tmesh), num_modes(model), num_structures(model))
    raman_noise = zeros(typeof(sqrt(ħ*ω0)), Base.size(raman_response))
    for (i, structure) in enumerate(model.waveguide.structures)
        for (j, mode) in enumerate(structure.medium.modes)
            if mode.has_raman
                tau1 = structure.medium.material.raman.tau1
                tau2 = structure.medium.material.raman.tau2
                raman_timeresponse = Complex.((tmesh .> 0ps) .* (tau1.^2 + tau2.^2)./tau1./tau2.^2 .* exp.(-tmesh/tau2).*sin.(tmesh/tau1)*oneunit(1.0ps))
                raman_response[:, j ,i] = ifft!(fftshift(raman_timeresponse))
                raman_noise[:, j, i] = sqrt.(2 * ħ * ω0 * abs.(imag.(view(raman_response, :, j, i))) .* (bose_distribution.(abs.(ω-ω0)+ω[2]-ω[1]) .+ (ω .< ω0)))
            end
        end
    end
    return raman_response, raman_noise
end
@inline bose_distribution(Ω) = one(Ω) / (exp(ħ_over_kBT * Ω) - 1)
