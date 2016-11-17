using StochasticDiffEq, DiffEqDevTools, DiffEqProblemLibrary
using Base.Test

const TEST_PLOT = false

#SDE
@time @testset "Linear SDE Tests" begin include("sde/sde_linear_tests.jl") end
@time @testset "Two-dimensional Linear SDE Tests" begin include("sde/sde_twodimlinear_tests.jl") end
@time @testset "Additive SDE Tests" begin include("sde/sde_additive_tests.jl") end
@time @testset "Rossler Order Tests" begin include("sde/sde_rosslerorder_tests.jl") end
@time @testset "SDE Convergence Tests" begin include("sde/sde_convergence_tests.jl") end
@time @testset "SDE Number Type Tests" begin include("sde/sde_numbertype_tests.jl") end
@time @testset "Oval2" begin include("oval2_test.jl") end

#Adaptive SDE
#=
println("Adaptive SDE Linear Tests")
@time @test include("sde/sde_linearadaptive_tests.jl")
println("Adaptive SDE Distribution Test")
@time @test include("sde/sde_adaptivedistribution_tests.jl")
println("Multiple Dimension Linear Adaptive Test")
@time @test include("sde/sde_twodimlinearadaptive_tests.jl")
println("SDE Autostepsize Test")
@time @test include("sde/sde_autostepsize_test.jl")
println("SDE Additive Lorenz Attractor Test")
@time @test include("sde/sde_lorenzattractor_tests.jl")
=#
