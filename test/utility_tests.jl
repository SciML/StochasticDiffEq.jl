using StochasticDiffEq, LinearAlgebra, SparseArrays, Random, Test, DiffEqOperators
using StochasticDiffEq: WOperator, set_gamma!, calc_W!

@testset "Derivative Utilities" begin
  @testset "WOperator" begin
    Random.seed!(123)
    y = zeros(2); b = rand(2)
    mm = rand(2, 2)
    for _J in [rand(2, 2), Diagonal(rand(2))]
      _Ws = [-mm + 2.0 * _J, -mm/2.0 + _J]
      for inplace in (true, false), (_W, W_transform) in zip(_Ws, [false, true])
        W = WOperator(mm, 1.0, DiffEqArrayOperator(_J), inplace, transform=W_transform)
        set_gamma!(W, 2.0)
        @test convert(AbstractMatrix,W) ≈ _W
        @test W * b ≈ _W * b
        mul!(y, W, b); @test y ≈ _W * b
      end
    end
  end

  @testset "calc_W!" begin
    A = [-1.0 0.0; 0.0 -0.5]; σ = [0.9 0.0; 0.0 0.8]
    mm = [2.0 0.0; 0.0 1.0]
    u0 = [1.0, 1.0]; tmp = zeros(2)
    tspan = (0.0,1.0); dt = 0.01
    concrete_W = -mm + dt * A

    # Out-of-place
    _f = (u,p,t) -> A*u; _g = (u,p,t) -> σ*u
    fun = SDEFunction(_f, _g;
                      mass_matrix=mm,
                      jac=(u,p,t) -> A)
    prob = SDEProblem(fun, _g, u0, tspan)
    integrator = init(prob, ImplicitEM(theta=1); adaptive=false, dt=dt)
    J, W = calc_W!(integrator.cache.nlsolver, integrator, integrator.cache, dt, false)
    @test convert(AbstractMatrix, W) ≈ concrete_W
    @test W \ u0 ≈ concrete_W \ u0

    # In-place
    _f = (du,u,p,t) -> mul!(du,A,u); _g = (du,u,p,t) -> mul!(du,σ,u)
    fun = SDEFunction(_f, _g;
                      mass_matrix=mm,
                      jac_prototype=DiffEqArrayOperator(A))
    prob = SDEProblem(fun, _g, u0, tspan)
    integrator = init(prob, ImplicitEM(theta=1); adaptive=false, dt=dt)
    calc_W!(integrator.cache.nlsolver, integrator, integrator.cache, dt, false)
    @test convert(AbstractMatrix, integrator.cache.nlsolver.cache.W) ≈ concrete_W
    ldiv!(tmp, lu!(integrator.cache.nlsolver.cache.W), u0); @test tmp ≈ concrete_W \ u0
  end

  @testset "Implicit solver with lazy W" begin
    A = sparse([-1.0 0.0; 0.0 -0.5]); σ = sparse([0.9 0.0; 0.0 0.8])
    mm = sparse([2.0 0.0; 0.0 1.0])
    u0 = [1.0, 1.0]; tspan = (0.0,1.0)

    _f = (u,p,t) -> t*(A*u); _f_ip = (du,u,p,t) -> lmul!(t,mul!(du,A,u))
    _g = (u,p,t) -> σ*u; _g_ip = (du,u,p,t) -> mul!(du,σ,u)
    prob1 = SDEProblem(SDEFunction(_f, _g; mass_matrix=mm), _g, u0, tspan)
    prob2 = SDEProblem(SDEFunction(_f, _g; mass_matrix=mm, jac=(u,p,t) -> t*A), _g, u0, tspan)
    prob1_ip = SDEProblem(SDEFunction(_f_ip, _g_ip; mass_matrix=mm), _g_ip, u0, tspan)
    jac_prototype=DiffEqArrayOperator(similar(A); update_func=(J,u,p,t) -> (J .= t .* A; J))
    prob2_ip = SDEProblem(SDEFunction(_f_ip, _g_ip; mass_matrix=mm, jac_prototype=jac_prototype), _g_ip, u0, tspan)

    for Alg in [ImplicitEM, ISSEM]
      println(Alg)
      Random.seed!(0); sol1 = solve(prob1, Alg(theta=1); adaptive=false, dt=0.01)
      Random.seed!(0); sol2 = solve(prob2, Alg(theta=1); adaptive=false, dt=0.01)
      @test sol1(1.0) ≈ sol2(1.0) rtol=1e-4
      Random.seed!(0); sol1_ip = solve(prob1_ip, Alg(theta=1); adaptive=false, dt=0.01)
      Random.seed!(0); sol2_ip = solve(prob2_ip, Alg(theta=1,linsolve=LinSolveFactorize(lu)); adaptive=false, dt=0.01)
      @test sol1_ip(1.0) ≈ sol2_ip(1.0) rtol=1e-4
    end

    σ = 1.0
    _g = (u,p,t) -> σ; _g_ip = (du,u,p,t) -> (du .= σ)
    prob1 = SDEProblem(SDEFunction(_f, _g), _g, u0, tspan)
    prob2 = SDEProblem(SDEFunction(_f, _g; jac=(u,p,t) -> t*A), _g, u0, tspan)
    prob1_ip = SDEProblem(SDEFunction(_f_ip, _g_ip), _g_ip, u0, tspan)
    prob2_ip = SDEProblem(SDEFunction(_f_ip, _g_ip; jac_prototype=jac_prototype), _g_ip, u0, tspan)

    println(SKenCarp)
    Random.seed!(0); sol1 = solve(prob1, SKenCarp(); adaptive=false, dt=0.01)
    Random.seed!(0); sol2 = solve(prob2, SKenCarp(); adaptive=false, dt=0.01)
    @test sol1(1.0) ≈ sol2(1.0) rtol=1e-4
    Random.seed!(0); sol1_ip = solve(prob1_ip, SKenCarp(); adaptive=false, dt=0.01)
    Random.seed!(0); sol2_ip = solve(prob2_ip, SKenCarp(linsolve=LinSolveFactorize(lu)); adaptive=false, dt=0.01)
    @test sol1_ip(1.0) ≈ sol2_ip(1.0) rtol=1e-3
  end
end
