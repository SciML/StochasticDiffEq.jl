using SafeTestsets

const LONGER_TESTS = false

const GROUP = get(ENV, "GROUP", "All")

const is_APPVEYOR = Sys.iswindows() && haskey(ENV,"APPVEYOR")

@time begin
  if GROUP == "All" || GROUP == "Interface"
    @time @safetestset "First Rand Tests" begin include("first_rand_test.jl") end
    @time @safetestset "Inference Tests" begin include("inference_test.jl") end
    @time @safetestset "Linear RODE Tests" begin include("rode_linear_tests.jl") end
    @time @safetestset "Complex Number Tests" begin include("complex_tests.jl") end
    @time @safetestset "Static Array Tests" begin include("static_array_tests.jl") end
    @time @safetestset "Noise Type Tests" begin include("noise_type_test.jl") end
    @time @safetestset "Mass matrix tests" begin include("mass_matrix_tests.jl") end
    #@time @safetestset "Sparse Jacobian tests" begin include("sparsediff_tests.jl") end
    @time @safetestset "Outofplace Arrays Tests" begin include("outofplace_arrays.jl") end
    @time @safetestset "tdir Tests" begin include("tdir_tests.jl") end
    @time @safetestset "tstops Tests" begin include("tstops_tests.jl") end
    @time @safetestset "saveat Tests" begin include("saveat_tests.jl") end
    @time @safetestset "Oval2" begin include("oval2_test.jl") end
  end

  if GROUP == "All" || GROUP == "Interface2"
    @time @safetestset "Basic Tau Leaping Tests" begin include("tau_leaping.jl") end
    @time @safetestset "Linear SDE Tests" begin include("sde/sde_linear_tests.jl") end
    @time @safetestset "Two-dimensional Linear SDE Tests" begin include("sde/sde_twodimlinear_tests.jl") end
    @time @safetestset "Element-wise Tolerances Tests" begin include("tolerances_tests.jl") end
    @time @safetestset "Zero'd Noise Tests" begin include("zerod_noise_test.jl") end
    @time @safetestset "Scalar Tests" begin include("scalar_noise.jl") end
    @time @safetestset "Stiffness Detection Test" begin include("stiffness_detection_test.jl") end
    @time @safetestset "Adaptive SDE Linear Tests" begin include("adaptive/sde_linearadaptive_tests.jl") end
  end

  if GROUP == "All" || GROUP == "Interface3"
    @time @safetestset "Composite Tests" begin include("composite_algorithm_test.jl") end
    @time @safetestset "Events Tests" begin include("events_test.jl") end
    @time @safetestset "Cache Tests" begin include("cache_test.jl") end
    @time @safetestset "Adaptive Complex Mean Test" begin include("adaptive/sde_complex_adaptive_mean_test.jl") end
    @time @safetestset "Utility Tests" begin include("utility_tests.jl") end
    @time @safetestset "Non-diagonal SDE Tests" begin include("nondiagonal_tests.jl") end
    @time @safetestset "No Index Tests" begin include("noindex_tests.jl") end
    @time @safetestset "Multiple Dimension Linear Adaptive Test" begin include("adaptive/sde_twodimlinearadaptive_tests.jl") end
    @time @safetestset "Autostepsize Test" begin include("adaptive/sde_autostepsize_test.jl") end
    @time @safetestset "Additive Lorenz Attractor Test" begin include("adaptive/sde_lorenzattractor_tests.jl") end
  end

  if !is_APPVEYOR && (GROUP == "All" || GROUP == "AlgConvergence")
    @time @safetestset "Convergence Tests" begin include("sde/sde_convergence_tests.jl") end
  end

  if !is_APPVEYOR && (GROUP == "All" || GROUP == "AlgConvergence2")
    @time @safetestset "IIF Convergence Tests" begin include("iif_methods.jl") end
    @time @safetestset "Cummutative Noise Methods Tests" begin include("commutative_tests.jl") end
    @time @safetestset "Multivariate Geometric Tests" begin include("multivariate_geometric.jl") end
  end

  if !is_APPVEYOR && (GROUP == "All" || GROUP == "AlgConvergence3")
    @time @safetestset "Rossler Order Tests" begin include("sde/sde_rosslerorder_tests.jl") end
    @time @safetestset "ODE Convergence Regression Tests" begin include("ode_convergence_regression.jl") end
    @time @safetestset "Additive SDE Tests" begin include("sde/sde_additive_tests.jl") end
    @time @safetestset "Split Tests" begin include("split_tests.jl") end
    @time @safetestset "Stratonovich Convergence Tests" begin include("stratonovich_convergence_tests.jl") end
  end

  if !is_APPVEYOR && (GROUP == "All" || GROUP == "WeakConvergence")
    @time @safetestset "Roessler weak SRK Tests" begin include("weak_convergence/srk_weak_final.jl") end
    @time @safetestset "OOP Weak Convergence Tests" begin include("weak_convergence/oop_weak.jl") end
    @time @safetestset "IIP Weak Convergence Tests" begin include("weak_convergence/iip_weak.jl") end
    @time @safetestset "Additive Weak Convergence Tests" begin include("weak_convergence/additive_weak.jl") end
  end
end
