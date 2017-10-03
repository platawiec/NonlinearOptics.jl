using Base.Test
using NonlinearOptics

@time @testset "NLSE Tests" begin include("nlse/nlse_tests.jl") end
