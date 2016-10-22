using StochasticDiffEq
using Base.Test

const TEST_PLOT = false

# write your own tests here
#SDE
println("Linear SDE Tests")
@time @test include("sde/sde_linear_tests.jl")
println("Two-dimensional Linear SDE Tests")
@time @test include("sde/sde_twodimlinear_tests.jl")
println("Additive SDE Tests")
@time @test include("sde/sde_additive_tests.jl")
println("Rossler Order Tests")
@time @test include("sde/sde_rosslerorder_tests.jl")
println("SDE Convergence Tests")
@time @test include("sde/sde_convergence_tests.jl")
println("SDE Number Type Tests")
@time @test include("sde/sde_numbertype_tests.jl")
println("Oval2")
@time @test include("oval2_test.jl")

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
