using NonlinearOptics

# Define a soliton problem
prob = prob_bright_soliton

# Solve it
sol = solve(prob, SymmetrizedSplitStep())

prob = prob_dispersive_pulse
sol = solve(prob, SymmetrizedSplitStep())

prob = prob_LL
sol = solve(prob, SymmetrizedSplitStep())
