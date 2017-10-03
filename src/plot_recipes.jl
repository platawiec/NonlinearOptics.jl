@recipe function f(sol::NLSESolution;plot_analytic=false)
    if sol.tslocation == 0
        @series begin
            seriestype := :path
            sol.prob.zmesh, abs2.(sol.u[end])
        end
    else
        out = Any[]
        @series begin
            seriestype --> :surface
            layout --> length(out)
            sol.prob.zmesh, sol.prob.timeseries, sol.u
        end
    end
end
