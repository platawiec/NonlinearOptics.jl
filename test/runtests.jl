using Base.Test
using NonlinearOptics
import NonlinearOptics

detected_tests = filter(name->startswith(name, "test_") && endswith(name, ".jl"),
                        readdir(@__DIR__))

@time @testset "NonlinearOptics" begin
    @time @testset "$(test[6:(end-3)])" for test in detected_tests
        include(test)
    end
end
