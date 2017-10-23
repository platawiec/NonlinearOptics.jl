@recipe function f(sol::DynamicNLSESolution;plot_analytic=false)
    ω = sol.prob.ω
    zend = sol.prob.prob.tspan[2]
    spectrum = abs2.(FT(sol, zend))
    layout := (2, 1)
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
