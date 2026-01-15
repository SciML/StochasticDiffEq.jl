using StochasticDiffEq, OrdinaryDiffEqCore
using StochasticDiffEq: SDEVerbosity
using SciMLLogging: SciMLLogging, AbstractMessageLevel
using Test

@testset "SDEVerbosity Construction" begin
    # Test default constructor
    @testset "Default constructor" begin
        v = SDEVerbosity()
        @test v.ode_verbosity isa OrdinaryDiffEqCore.ODEVerbosity
        @test v.noise_evaluation isa SciMLLogging.Silent
        @test v.adaptive_timestepping isa SciMLLogging.Silent
        @test v.dt_NaN isa SciMLLogging.WarnLevel
        @test v.function_NaN isa SciMLLogging.WarnLevel
    end

    # Test preset constructors
    @testset "Preset: None" begin
        v = SDEVerbosity(None())
        @test v.ode_verbosity isa OrdinaryDiffEqCore.ODEVerbosity
        @test all(getfield(v, f) isa SciMLLogging.Silent
        for f in fieldnames(SDEVerbosity) if f != :ode_verbosity)
    end

    @testset "Preset: Minimal" begin
        v = SDEVerbosity(Minimal())
        @test v.noise_evaluation isa SciMLLogging.Silent
        @test v.adaptive_timestepping isa SciMLLogging.Silent
        @test v.dt_NaN isa SciMLLogging.WarnLevel
        @test v.function_NaN isa SciMLLogging.WarnLevel
    end

    @testset "Preset: Standard" begin
        v = SDEVerbosity(Standard())
        @test v isa SDEVerbosity
        @test v.dt_NaN isa SciMLLogging.WarnLevel
        @test v.function_NaN isa SciMLLogging.WarnLevel
    end

    @testset "Preset: Detailed" begin
        v = SDEVerbosity(Detailed())
        @test v.noise_evaluation isa SciMLLogging.InfoLevel
        @test v.adaptive_timestepping isa SciMLLogging.InfoLevel
        @test v.dt_NaN isa SciMLLogging.WarnLevel
        @test v.function_NaN isa SciMLLogging.WarnLevel
    end

    @testset "Preset: All" begin
        v = SDEVerbosity(All())
        @test v.noise_evaluation isa SciMLLogging.InfoLevel
        @test v.adaptive_timestepping isa SciMLLogging.InfoLevel
        @test v.dt_NaN isa SciMLLogging.WarnLevel
        @test v.function_NaN isa SciMLLogging.WarnLevel
    end

    # Test group-level constructor
    @testset "Group-level construction" begin
        v = SDEVerbosity(sde_specific = InfoLevel())
        @test v.noise_evaluation isa SciMLLogging.InfoLevel
        @test v.adaptive_timestepping isa SciMLLogging.InfoLevel
        @test v.dt_NaN isa SciMLLogging.InfoLevel
        @test v.function_NaN isa SciMLLogging.InfoLevel
    end

    # Test individual field construction
    @testset "Individual field construction" begin
        v = SDEVerbosity(
            noise_evaluation = InfoLevel(),
            dt_NaN = WarnLevel(),
            function_NaN = DebugLevel()
        )
        @test v.noise_evaluation isa SciMLLogging.InfoLevel
        @test v.dt_NaN isa SciMLLogging.WarnLevel
        @test v.function_NaN isa SciMLLogging.DebugLevel
        @test v.adaptive_timestepping isa SciMLLogging.Silent
    end

    # Test mixed group and individual settings
    @testset "Mixed group and individual" begin
        v = SDEVerbosity(
            sde_specific = InfoLevel(),
            dt_NaN = WarnLevel()  # Override
        )
        @test v.noise_evaluation isa SciMLLogging.InfoLevel
        @test v.adaptive_timestepping isa SciMLLogging.InfoLevel
        @test v.dt_NaN isa SciMLLogging.WarnLevel  # Overridden
    end

    # Test ODE verbosity passthrough
    @testset "ODE verbosity" begin
        ode_v = OrdinaryDiffEqCore.ODEVerbosity(Detailed())
        v = SDEVerbosity(ode_verbosity = ode_v)
        @test v.ode_verbosity === ode_v
    end

    # Test invalid arguments
    @testset "Invalid arguments" begin
        @test_throws ArgumentError SDEVerbosity(invalid_field = InfoLevel())
        @test_throws ArgumentError SDEVerbosity(sde_specific = "not a level")
        @test_throws ArgumentError SDEVerbosity(noise_evaluation = 42)
    end
end

@testset "SDEVerbosity Property Access" begin
    # Test direct field access (SDE-specific)
    @testset "Direct SDE field access" begin
        v = SDEVerbosity(noise_evaluation = InfoLevel())
        @test v.noise_evaluation isa SciMLLogging.InfoLevel
    end

    # Test nested ODE field access via getproperty
    @testset "Nested ODE field access" begin
        ode_v = OrdinaryDiffEqCore.ODEVerbosity(dt_NaN = WarnLevel())
        v = SDEVerbosity(ode_verbosity = ode_v)
        @test v.dt_NaN isa SciMLLogging.AbstractMessageLevel
        @test v.max_iters isa SciMLLogging.AbstractMessageLevel
    end

    # Test that invalid field access errors
    @testset "Invalid field access" begin
        v = SDEVerbosity()
        @test_throws ErrorException v.nonexistent_field
    end
end

@testset "SDEVerbosity in solve" begin
    # Define a simple SDE problem
    function f(du, u, p, t)
        du[1] = -u[1]
    end
    function g(du, u, p, t)
        du[1] = 0.1u[1]
    end
    prob = SDEProblem(f, g, [1.0], (0.0, 1.0))

    # Test Bool verbose conversion
    @testset "Bool verbose argument" begin
        sol = solve(prob, EM(); verbose = true, dt = 0.01)
        @test sol.retcode == ReturnCode.Success

        sol = solve(prob, EM(); verbose = false, dt = 0.01)
        @test sol.retcode == ReturnCode.Success
    end

    # Test AbstractVerbosityPreset conversion
    @testset "Preset verbose argument" begin
        sol = solve(prob, EM(); verbose = Standard(), dt = 0.01)
        @test sol.retcode == ReturnCode.Success

        sol = solve(prob, EM(); verbose = None(), dt = 0.01)
        @test sol.retcode == ReturnCode.Success
    end

    # Test SDEVerbosity directly
    @testset "SDEVerbosity argument" begin
        v = SDEVerbosity(noise_evaluation = InfoLevel())
        sol = solve(prob, EM(); verbose = v, dt = 0.01)
        @test sol.retcode == ReturnCode.Success
    end

    # Test that integrator has correct verbosity type
    @testset "Integrator verbosity type" begin
        integrator = init(prob, EM(); verbose = true, dt = 0.01)
        @test integrator.opts.verbose isa SDEVerbosity

        integrator = init(prob, EM(); verbose = Standard(), dt = 0.01)
        @test integrator.opts.verbose isa SDEVerbosity

        v = SDEVerbosity(dt_NaN = WarnLevel())
        integrator = init(prob, EM(); verbose = v, dt = 0.01)
        @test integrator.opts.verbose === v
    end
end
