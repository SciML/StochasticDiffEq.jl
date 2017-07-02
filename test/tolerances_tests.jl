using StochasticDiffEq, DiffEqProblemLibrary, Base.Test

#=
f = (t,u) -> begin
  du = similar(u)
  prob_sde_2Dlinear.f(t,u,du)
  return du
end
σ = (t,u) -> begin
  du = similar(u)
  prob_sde_2Dlinear.g(t,u,du)
  return du
end

test_probs = [prob_sde_2Dlinear,
	      SDEProblem(f,σ,prob_sde_2Dlinear.u0,prob_sde_2Dlinear.tspan)]
=#

test_probs = [prob_sde_2Dlinear]
adaptive_test_algs = [SRI(),SRIW1(),SRA1(),SRA()]
fixed_test_algs = [RKMil()]

for alg in adaptive_test_algs, prob in test_probs
  srand(100)
  sol = solve(prob,alg)

  # Vector of element-wise absolute tolerances
  srand(100)
  sol2 = solve(prob,alg;abstol=fill(1e-2,4,2))
  
  @test sol.t == sol2.t && sol.u == sol2.u

  # Vector of element-wise relative tolerances
  srand(100)
  sol2 = solve(prob,alg;reltol=fill(1e-2,4,2))

  @test sol.t == sol2.t && sol.u == sol2.u
end

for alg in fixed_test_algs, prob in test_probs
  srand(100)
  sol = solve(prob,alg;dt=0.1)

  # Vector of element-wise absolute tolerances
  srand(100)
  sol2 = solve(prob,alg;dt=0.1,abstol=fill(1e-2,4,2))
  
  @test sol.t == sol2.t && sol.u == sol2.u

  # Vector of element-wise relative tolerances
  srand(100)
  sol2 = solve(prob,alg;dt=0.1,reltol=fill(1e-2,4,2))

  @test sol.t == sol2.t && sol.u == sol2.u
end
