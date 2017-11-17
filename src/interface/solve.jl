function solve(model::AbstractModel,
               probtype::AbstractExperiment,
               alg::AbstractODEAlgorithm=DP5();
               kwargs...)
    prob = build_problem(model, probtype; kwargs...)
    sol = solve(prob, alg; kwargs...)
    return sol
end

function solve(prob::AbstractNLOProblem, alg::AbstractODEAlgorithm=DP5(); kwargs...)
    sol = solve(prob.prob, alg; kwargs...)
    DynamicNLOSolution(sol, prob)
end
function solve(prob::DynamicIkedaProblem, alg::AbstractODEAlgorithm=DP5(); kwargs...)
    sol = solve(prob.prob, alg; callback=prob.ikeda_callback, kwargs...)
    DynamicNLOSolution(sol, prob)
end
