# preferred units are THz, ps, and m
freq = collect(linspace(150, 250, 250))
neff_data = collect(linspace(3.0, 1.0, 250))
neff = OpticalAttr(freq, neff_data, "n_eff")
corefrac = OpticalAttr(1.0, "Core Fraction")
Aeff = OpticalAttr(1e-12, "Aeff (m²)")
loss = OpticalAttr(0.0, "Loss (1/m)")
coupling = OpticalAttr(0.6, "Coupling")

mode_TM = Mode(neff, Aeff, loss, corefrac, coupling, :TM, true, true)
mode_TE = Mode(neff, Aeff, loss, corefrac, coupling, :TE, true, true)

# 10 mm long waveguide in 100 orientation
wg_simple = Waveguide(1e-3, 100, Silicon)
add_mode!(wg_simple, mode_TM)
add_mode!(wg_simple, mode_TE)
add_interaction!(wg_simple, mode_TM, mode_TE)
# Note: pulse width is in 1/THz = fs
source_simple = PulsedLaser(Frequency(200.), 10.0, 0.020)
model_NLSE = Model(source_simple, [wg_simple])
# Laser answers: δ₀ can change during experiment, ω₀, and Ein? All can be
# parameters of cw laser. Solver ?s: How many dispersive orders
# should we use (keyword, default=3)? Initial conditions for NLSE? How do we
# determine initial conditions otherwise?
#solve(model_simple, DynamicLL()) (DynamicLL, DynamicIkeda, SSLL, SSIkeda, NLSE)

res_simple = CircularResonator(300e-6, SiO2)
add_mode!(res_simple, mode)
source_CW = CWLaser(Frequency(200), 0.05, 1.0)
model_LL = Model(source_CW, [res_simple, res_simple])

sol_NLSE = solve(model_NLSE, StochasticNLSE(); tpoints=2^12)

R = CubicRamanTensor()
ETensor = CubicElectronicTensor(0.1)
