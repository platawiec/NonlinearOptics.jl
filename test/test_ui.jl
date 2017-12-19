freq = collect(linspace(150, 250, 50))THz
neff_data = collect(linspace(3.0, 1.0, 50))
neff = OpticalAttr(freq, neff_data, "n_eff")
corefrac = OpticalAttr(1.0, "Core Fraction")
Aeff = OpticalAttr(1e-12m^2, "Aeff (m²)")
loss = OpticalAttr(0.0m^-1, "Loss (1/m)")
coupling = OpticalAttr(0.6, "Coupling")

mode_TM = Mode(neff, Aeff, loss, corefrac, coupling, :TM, true, true)
mode_TE = Mode(neff, Aeff, loss, corefrac, coupling, :TE, true, true)
mode_TM2 = Mode(neff, Aeff, loss, corefrac, coupling, :TM, true, true)

wg_simple = Waveguide(Silicon)
add_straight!(wg_simple, 1e-3m)
add_turn!(wg_simple, 90°, 1e-3m)
add_mode!(wg_simple, mode_TM)
add_mode!(wg_simple, mode_TE)
add_mode!(wg_simple, mode_TM2)
@test mode_TM.has_raman
@test !(mode_TM.has_interaction)
@test !(mode_TE.has_interaction)
@test !NonlinearOptics.interacts_with_other_mode(mode_TM)
@test NonlinearOptics.get_polarization(mode_TM) == 1
@test NonlinearOptics.get_polarization(mode_TE) == 2

add_interaction!(wg_simple, mode_TM, mode_TE)
@test mode_TM.has_interaction
@test mode_TE.has_interaction
@test NonlinearOptics.interacts_with_other_mode(mode_TM)
@test !NonlinearOptics.interacts_with_other_mode(mode_TM2)

# link waveguides
