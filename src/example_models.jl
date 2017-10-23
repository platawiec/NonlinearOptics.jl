# preferred units are THz
freq = collect(linspace(150, 250, 250))
neff_data = collect(linspace(3.0, 1.0, 250))
neff = OpticalAttr(freq, neff_data, "n_eff")
corefrac = OpticalAttr(1.0, "Core Fraction")
Aeff = OpticalAttr(1e-12, "Aeff (m²)")
loss = OpticalAttr(1e-3, "Loss (1/m)")
coupling = OpticalAttr(0.6, "Coupling")

mode = Mode(neff, Aeff, loss, corefrac, coupling)

# 10 mm long waveguide in 100 orientation
wg_simple = Waveguide(1e-3, 100, SiO2)
add_mode!(wg_simple, mode)
# Note: pulse width is in 1/THz = fs
source_simple = PulsedLaser(Frequency(200), 1.0, 1.0, 0.020)
model_NLSE = Model(source_simple, wg_simple)
# Laser answers: δ₀ can change during experiment, ω₀, and Ein? All can be
# parameters of cw laser. Solver ?s: How many dispersive orders
# should we use (keyword, default=3)? Initial conditions for NLSE? How do we
# determine initial conditions otherwise?
#solve(model_simple, DynamicLL()) (DynamicLL, DynamicIkeda, SSLL, SSIkeda, NLSE)

res_simple = CircularResonator(300e-6, SiO2)
add_mode!(res_simple, mode)
source_CW = CWLaser(Frequency(200), 0.05, 1.0)
model_LL = Model(source_CW, res_simple)

# Second-order Soliton model from J. Lightwave Tech. Vol. 25 No. 12 Dec. 2007
# Johan Hult
model_ToyNLSE = ToyModel(;betacoeff=[0.0, 0.0, -0.01], nonlinearcoeff=0.01,
                         power_in=1.24e3, pulsetime=0.1, length=2.0)
