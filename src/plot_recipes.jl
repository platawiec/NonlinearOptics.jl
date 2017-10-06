@recipe function f(sol::NLSESolution;plot_analytic=false)
    if sol.tslocation == 0
        layout := @layout [intime
                           inspectrum]
        @series begin
            subplot := 1
            seriestype := :path
            sol.prob.zmesh, abs2.(sol.u[end])
        end
        dω = 2pi/maximum(prob.zmesh)
        ω_max = dktilde/2.0*(length(prob.zmesh)-1)
        ω = collect(-ktilde_max:dktilde:ktilde_max)
        spectrum = abs2.(fftshift(fft(sol.u[end])))
        spectrum_dB = 10*log10.(spectrum/maximum(spectrum))
        @series begin
            subplot := 2
            seriestype := :path
            ω, spectrum_dB
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
