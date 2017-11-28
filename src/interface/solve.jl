function solve(model::AbstractModel,
               probtype::AbstractDynamicExperiment,
               alg::AbstractODEAlgorithm=DP5();
               kwargs...)
    prob = build_problem(model, probtype; kwargs...)
    sol = solve(prob, alg; kwargs...)
    return sol
end

function solve(model::AbstractModel,
               probtype::AbstractStochasticExperiment,
               alg::AbstractSDEAlgorithm=SRIW1();
               kwargs...)
    prob = build_problem(model, probtype; kwargs...)
    sol = solve(prob, alg; kwargs...)
    return sol
end

function solve(prob::AbstractNLOProblem, alg; kwargs...)
    sol = DiffEqBase.solve(prob.prob, alg; kwargs...)
    return DynamicNLOSolution(sol, prob)
end
function solve(prob::DynamicIkedaProblem, alg; kwargs...)
    sol = DiffEqBase.solve(prob.prob, alg; callback=prob.ikeda_callback, kwargs...)
    DynamicNLOSolution(sol, prob)
end
