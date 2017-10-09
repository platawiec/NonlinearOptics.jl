using Base.Test
using NonlinearOptics

@testset "NonlinearOptics" begin
    for (root, dirs, files) in walkdir(@__DIR__)
        @time @testset "$dir" for dir in dirs
            @time @testset "$file" for file in files include(file) end
        end
    end
end



@time @testset "NLSE Tests" begin include("nlse/nlse_tests.jl") end
