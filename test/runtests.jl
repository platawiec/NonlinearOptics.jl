using NonlinearOptics
using Base.Test

println("NLSE Tests")
@time @test include("nlse/nlse_tests.jl")
