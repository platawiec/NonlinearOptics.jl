@recipe function f(sol::NLSESolution;plot_analytic=false)
    if sol.tslocation == 0
        layout := (2, 1)
        @series begin
            subplot := 1
            seriestype := :path
            sol.prob.zmesh, abs2.(sol[end])
        end
        dω = 2pi/maximum(sol.prob.zmesh)
        ω_max = dω/2.0*(length(sol.prob.zmesh)-1)
        ω = collect(-ω_max:dω:ω_max)
        spectrum = abs2.(fftshift(fft(sol[end])))
        spectrum_dB = 10*log10.(spectrum/maximum(spectrum))
        @series begin
            subplot := 2
            seriestype := :path
            fill := (:auto, minimum(spectrum_dB))
            ω, spectrum_dB
        end
    else
        out = Any[]
        @series begin
            ylabel --> ("(dB)")
            title --> "Spectrum"
            seriestype --> :surface
            layout --> length(out)
            sol.prob.zmesh, sol.prob.timeseries, sol.u
        end
    end
end
