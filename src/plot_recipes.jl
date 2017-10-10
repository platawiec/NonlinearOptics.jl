@recipe function f(sol::NLSESolution;plot_analytic=false, FROG_plot=false)
    dω = 2pi/maximum(sol.prob.zmesh)
    ω_max = dω/2.0*(length(sol.prob.zmesh)-1)
    ω = collect(-ω_max:dω:ω_max)
    spectrum_dB = []
    for u in sol.u
        spectrum = abs2.(fftshift(fft(u)))
        push!(spectrum_dB, 10*log10.(spectrum/maximum(spectrum)))
    end
    if !FROG_plot
        layout := (2, 1)
        @series begin
            title := "Intensity"
            subplot := 1
            seriestype := :path
            sol.prob.zmesh, abs2.(sol[end])
        end
        @series begin
            ylabel := ("(dB)")
            title := "Spectrum"
            subplot := 2
            seriestype := :path
            fill := (:auto, minimum(spectrum_dB[end]))
            ω, spectrum_dB[end]
        end
    else
        # heatmap works for plotlyjs, but contour for gr
        @series begin
            seriestype := :heatmap
            fill := (true, :plasma)
            sol.t, ω, hcat(spectrum_dB...)
        end
    end
end

@recipe function f(optical_property::AbstractOpticalProperty)
    fs = frequency(optical_property)
    λs = wavelength(optical_property)
    @series begin
        seriestype := :path
        λs, get_property(optical_property)
    end
end
