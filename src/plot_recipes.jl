@recipe function f(sol::AbstractNLOSolution;plot_analytic=false, evolution=false)
    ω = sol.prob.ω
    zend = sol.prob.prob.tspan[2]
    spectrum = abs2.(FT(sol, zend))
    layout := (2, 1)
    if evolution
        zmesh = linspace(0, zend, 100)
        intensity = 10*log10.(1e-30+abs2.(sol.(zmesh)))
        intensity = hcat(intensity...)
        spectrum = 10*log10.(1e-30+abs2.(FT.(sol, zmesh)))
        spectrum = hcat(spectrum...)
        @series begin
            title := "Intensity"
            subplot := 1
            #ylims := (minimum(sol.prob.tmesh)/25, maximum(sol.prob.tmesh)/25)
            clims := (maximum(intensity)-80, maximum(intensity))
            zmesh, sol.prob.tmesh, intensity
        end
        @series begin
            title := "Spectrum"
            subplot := 2
            clims := (maximum(spectrum)-80, maximum(spectrum))
            zmesh, ω, spectrum
        end
    else
        @series begin
            title := "Intensity"
            subplot := 1
            seriestype := :path
            sol.prob.tmesh, abs2.(sol(zend))
        end
        @series begin
            ylabel := ("(dB)")
            title := "Spectrum"
            subplot := 2
            seriestype := :path
            ω, spectrum/maximum(spectrum)
        end
    end
end

@recipe function f(optical_property::AbstractOpticalAttr)
    fs = frequency(optical_property)/1e12 #plot in THz
    λs = wavelength(optical_property)*1e9 #plot in nm
    @series begin
        xlabel := "Frequency (THz)"
        ylabel := get_label(optical_property)
        seriestype := :path
        xmirror := true
        fs, get_property(optical_property)
    end
end
