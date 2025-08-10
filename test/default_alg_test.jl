using StochasticDiffEq, DiffEqNoiseProcess, Random, LinearAlgebra
using Test

Random.seed!(100)

# Simple SDE problem for testing
f(u, p, t) = 1.01 * u
g(u, p, t) = 0.87 * u
u0 = 1.0
tspan = (0.0, 1.0)

@testset "Default SDE Algorithm Tests" begin
    @testset "Basic DefaultSDEAlgorithm" begin
        prob = SDEProblem(f, g, u0, tspan)

        # Test with explicit DefaultSDEAlgorithm
        alg = DefaultSDEAlgorithm()
        sol = solve(prob, alg)
        @test sol.retcode == ReturnCode.Success
        @test length(sol.t) > 2

        # Test with nothing (should use default)
        sol2 = solve(prob, nothing)
        @test sol2.retcode == ReturnCode.Success
        @test length(sol2.t) > 2

        # Test solve without algorithm (should use default)
        sol3 = solve(prob)
        @test sol3.retcode == ReturnCode.Success
        @test length(sol3.t) > 2
    end

    @testset "DefaultAdaptiveSDEAlgorithm with hints" begin
        prob = SDEProblem(f, g, u0, tspan)

        # Test with additive noise hint
        alg = DefaultAdaptiveSDEAlgorithm(alg_hints = [:additive])
        sol = solve(prob, alg)
        @test sol.retcode == ReturnCode.Success

        # Test with stiff hint
        alg = DefaultAdaptiveSDEAlgorithm(alg_hints = [:stiff])
        sol = solve(prob, alg)
        @test sol.retcode == ReturnCode.Success

        # Test with commutative hint
        alg = DefaultAdaptiveSDEAlgorithm(alg_hints = [:commutative])
        sol = solve(prob, alg)
        @test sol.retcode == ReturnCode.Success
    end

    @testset "Additive noise problems" begin
        # Additive noise SDE
        f_add(u, p, t) = -0.5 * u
        g_add(u, p, t) = 0.1
        prob_add = SDEProblem(f_add, g_add, u0, tspan)

        sol = solve(prob_add, nothing)
        @test sol.retcode == ReturnCode.Success
        @test length(sol.t) > 2
    end

    @testset "System of SDEs" begin
        # System of SDEs
        f_sys(du, u, p, t) = begin
            du[1] = 1.01 * u[1]
            du[2] = -0.5 * u[2]
        end
        g_sys(du, u, p, t) = begin
            du[1] = 0.3 * u[1]
            du[2] = 0.2 * u[2]
        end
        u0_sys = [1.0, 2.0]
        prob_sys = SDEProblem(f_sys, g_sys, u0_sys, tspan)

        sol = solve(prob_sys, nothing)
        @test sol.retcode == ReturnCode.Success
        @test length(sol.t) > 2
        @test size(sol.u[end]) == size(u0_sys)
    end

    @testset "Stiff SDE problems" begin
        # Stiff problem
        f_stiff(u, p, t) = -100.0 * u
        g_stiff(u, p, t) = 0.1 * u
        prob_stiff = SDEProblem(f_stiff, g_stiff, u0, (0.0, 0.1))

        alg = DefaultAdaptiveSDEAlgorithm(alg_hints = [:stiff])
        sol = solve(prob_stiff, alg)
        @test sol.retcode == ReturnCode.Success
    end

    @testset "Algorithm selection logic" begin
        prob = SDEProblem(f, g, u0, tspan)

        # Test default_sde_alg_choice
        hints = default_sde_alg_choice(prob)
        @test isa(hints, Vector{Symbol})

        # Test with mass matrix (should add :stiff hint)
        function f_mass(du, u, p, t)
            du[1] = 1.01 * u[1]
            du[2] = -0.5 * u[2]
        end
        function g_mass(du, u, p, t)
            du[1] = 0.3 * u[1]
            du[2] = 0.2 * u[2]
        end
        M = [1.0 0.0; 0.0 2.0]
        f_mass_func = ODEFunction(f_mass, mass_matrix = M)
        prob_mass = SDEProblem(f_mass_func, g_mass, [1.0, 2.0], tspan)

        hints_mass = default_sde_alg_choice(prob_mass)
        @test :stiff in hints_mass
    end

    @testset "isdefaultalg check" begin
        alg1 = DefaultSDEAlgorithm()
        @test StochasticDiffEq.isdefaultalg(alg1) == true

        alg2 = DefaultAdaptiveSDEAlgorithm(alg_hints = [:additive])
        @test StochasticDiffEq.isdefaultalg(alg2) == true

        # Non-default algorithm
        alg3 = EM()
        @test StochasticDiffEq.isdefaultalg(alg3) == false
    end
end
