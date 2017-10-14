solve(model::Model, exp::DynamicLL)
end

"""
    build_problem()
Sets up a Lugatio-Lefever problem
"""
function build_problem(laser::CWLaser, res::AbstractResonator; dispersion_order=2)
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
