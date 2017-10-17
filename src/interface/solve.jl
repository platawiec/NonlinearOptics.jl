function solve(model::Model, exp::DynamicLL)
    prob = build_problem(model.laser, model.structure)

    solve(prob, alg=NLSE; kwargs...)
end

"""
    build_problem()
Sets up a Lugatio-Lefever problem
"""
function build_problem(laser::CWLaser, res::AbstractResonator;
                       dispersion_order=2, additional_terms=[])
    # TODO: support mode coupling
    # TODO: currently only supports one mode
    mode = res.modes[1]

    FSR = get_FSR(res, mode, laser.frequency)
    α = mode.linearloss
    γnl = get_nonlinearcoeff(res, mode, laser.frequency)
    L = circumference(res)
    Ein = sqrt(laser.power)
    detuning = laser.detuning
    beta = get_beta(mode, laser.frequency, dispersion_order)

    tmesh = linspace(-0.5, 0.5, 2^10)
    """Lugatio-Lefever evolution of cavity mean-field. From
    "Frequencu Comb Generation beyond the Lugatio-Lefever equation: multi-stability
    and super cavity solitons" T. Hansson, S. Wabnitz, Arxiv 1503.03274
    Parameters from Fig. 7a
    """
    prob_LL = NLSEProblem(N_LL, D_LL, u0_LL, τspan, tmesh)
end

function build_problem(laser::PulsedLaser, wg::Waveguide;
                       dispersion_order=2, additional_terms=[])
    # TODO: support mode coupling
    # TODO: currently only supports one mode
    wg = wg.modes[1]

    α = mode.linearloss
    γnl = get_nonlinearcoeff(res, mode, laser.frequency)
    L = wg.length

    beta = get_beta(mode, laser.frequency, dispersion_order)
    beta_coeff = beta ./ [factorial(n) for n=0:(length(beta)-1)]
    tmesh = linspace(-0.5, 0.5, 2^10)

    u0 = derive_pulse(laser)
    self_steepening(z, u) = 0.0
    raman_response(z, u) = 0.0
    if :self_steepening in additional_terms
        ω0 = getω(laser.frequency)
        self_steepening(z, u) = - γnl/ω0 * diff_cyclic(abs2(u)) / dt * u
    end
    if :raman_response in additional_terms
        raman_response(z, u) - 1im*γnl*raman_response(wg.material)
    end

    function N(z, u)
        (-α/2 + 1im * γnl * abs2(u)) * u + self_steepening(z, u) + raman_response(z, u)
    end
    function D(ω, u)
        1im * Poly(beta_coeff, :ω)
    end
    prob = NLSEProblem(N, D, u0, zspan, tmesh)
end

function diff_cyclic(A)
    B = similar(A)
    @inbounds B[1] = A[1]-A[end]
    for i=2:length(B)
        @inbounds B[i] = A[i]-A[i-1]
    end
    B
end

function diff_cyclic!(A)
    @inbounds tmp = A[1]-A[end]
    for i=1:(length(A)-1)
        @inbounds A[i] = A[i+1]-A[i]
    end
    @inbounds A[end] = tmp
    A
end
