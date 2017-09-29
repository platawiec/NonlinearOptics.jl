using NonlinearOptics

# Define a soliton problem
prob = prob_bright_soliton

# Solve it
println("Split-Step Fourier Method")
sol = solve(prob, NLSESplitStepFourier())
