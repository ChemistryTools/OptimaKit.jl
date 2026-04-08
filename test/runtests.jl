using Test
using OptimaJL
import OptimaJL: solve, solve!   # not exported; import explicitly for tests
using LinearAlgebra
using ForwardDiff

@testset "OptimaJL" begin
    include("test_canonicalizer.jl")
    include("test_newton.jl")
    include("test_solver.jl")
    include("test_sensitivity.jl")
    include("test_ad.jl")
end
