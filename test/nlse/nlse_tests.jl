using NonlinearOptics

# Define a soliton problem
prob = prob_bright_soliton
"""Example of a bright soliton problem construction
"""
# Solve it
sol = solve(prob, SymmetrizedSplitStep())
