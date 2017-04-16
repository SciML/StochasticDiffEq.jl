using StochasticDiffEq, DiffEqDevTools, DiffEqProblemLibrary
using Base.Test

tic()
const TEST_PLOT = false
const LONGER_TESTS = false
#SDE
@time @testset "First Rand Tests" begin include("sde/first_rand_test.jl") end
@time @testset "Linear SDE Tests" begin include("sde/sde_linear_tests.jl") end
@time @testset "Linear RODE Tests" begin include("rode_linear_tests.jl") end
@time @testset "Two-dimensional Linear SDE Tests" begin include("sde/sde_twodimlinear_tests.jl") end
@time @testset "Additive SDE Tests" begin include("sde/sde_additive_tests.jl") end
@time @testset "Rossler Order Tests" begin include("sde/sde_rosslerorder_tests.jl") end
@time @testset "Convergence Tests" begin include("sde/sde_convergence_tests.jl") end
@time @testset "Stratonovich Convergence Tests" begin include("stratonovich_convergence_tests.jl") end
LONGER_TESTS && @time @testset "Weak Convergence Tests" begin include("weak_convergence.jl") end
@time @testset "Number Type Tests" begin include("sde/sde_numbertype_tests.jl") end
@time @testset "Noise Type Tests" begin include("noise_type_test.jl") end
@time @testset "tdir Tests" begin include("tdir_tests.jl") end
@time @testset "tstops Tests" begin include("tstops_tests.jl") end
@time @testset "saveat Tests" begin include("saveat_tests.jl") end
@time @testset "Oval2" begin include("oval2_test.jl") end
@time @testset "Composite Tests" begin include("composite_algorithm_test.jl") end
@time @testset "Events Tests" begin include("events_test.jl") end
@time @testset "Cache Tests" begin include("cache_test.jl") end

#Adaptive SDE
@time @testset "Adaptive SDE Linear Tests" begin include("adaptive/sde_linearadaptive_tests.jl") end
@time @testset "Adaptive SDE Distribution Test" begin include("adaptive/sde_adaptivedistribution_tests.jl") end
@time @testset "Multiple Dimension Linear Adaptive Test" begin include("adaptive/sde_twodimlinearadaptive_tests.jl") end
@time @testset "Autostepsize Test" begin include("adaptive/sde_autostepsize_test.jl") end
@time @testset "Additive Lorenz Attractor Test" begin include("adaptive/sde_lorenzattractor_tests.jl") end
toc()
