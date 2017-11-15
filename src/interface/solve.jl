function solve(prob::AbstractNLOProblem, alg::AbstractODEAlgorithm=DP5(); kwargs...)
    sol = solve(prob.prob, alg; kwargs...)
    DynamicNLOSolution(sol, prob)
end
function solve(prob::DynamicIkedaProblem, alg::AbstractODEAlgorithm=DP5(); kwargs...)
    sol = solve(prob.prob, alg; callback=prob.ikeda_callback, kwargs...)
    DynamicNLOSolution(sol, prob)
end
