# preferred units are THz, ps, and m
freq = collect(linspace(150, 250, 250))THz
neff_data = collect(linspace(3.0, 1.0, 250))
neff = OpticalAttr(freq, neff_data, "n_eff")
corefrac = OpticalAttr(1.0, "Core Fraction")
Aeff = OpticalAttr(1e-12m^2, "Aeff (m²)")
loss = OpticalAttr(0.0m^-1, "Loss (1/m)")
coupling = OpticalAttr(0.6, "Coupling")

mode_TM = Mode(neff, Aeff, loss, corefrac, coupling, :TM, true, true)
mode_TE = Mode(neff, Aeff, loss, corefrac, coupling, :TE, true, true)

# 10 mm long waveguide in default (100) orientation
wg_simple = Waveguide(Silicon)
add_straight!(wg_simple, 1e-3m)
add_turn!(wg_simple, 90°, 1e-3m)
add_mode!(wg_simple, mode_TM)
add_mode!(wg_simple, mode_TE)
add_interaction!(wg_simple, mode_TM, mode_TE)
#add_medium!(wg_simple, Diamond)
#
# Note: pulse width is in 1/THz = fs
source_simple = PulsedLaser(200.0THz, 1.0W, 0.20ps)
model_NLSE = Model(source_simple, wg_simple)
# Laser answers: δ₀ can change during experiment, ω₀, and Ein? All can be
# parameters of cw laser. Solver ?s: How many dispersive orders
# should we use (keyword, default=3)? Initial conditions for NLSE? How do we
# determine initial conditions otherwise?
#solve(model_simple, DynamicLL()) (DynamicLL, DynamicIkeda, SSLL, SSIkeda, NLSE)

res_simple = CircularResonator(300e-6, SiO2)
add_mode!(res_simple, mode_TE)
source_CW = CWLaser(200THz, 0.05, 1.0)
model_LL = Model(source_CW, [res_simple, res_simple])

sol_NLSE = solve(model_NLSE, StochasticNLSE(); time_pts=2^5)
