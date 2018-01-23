using StochasticDiffEq, DiffEqProblemLibrary, Base.Test

#=
function f(u,p,t)
  du = similar(u)
  prob_sde_2Dlinear.f(du,u,p,t)
  return du
end
function σ(u,p,t)
  du = similar(u)
  prob_sde_2Dlinear.g(du,u,p,t)
  return du
end

probs = [prob_sde_2Dlinear,
         SDEProblem(f,σ,prob_sde_2Dlinear.u0,prob_sde_2Dlinear.tspan)]
=#

probs = [prob_sde_2Dlinear]
algs = [SRI(),SRIW1(),SRA1(),SRA(),RKMil()]

for alg in algs, prob in probs
  dt = typeof(alg)<:StochasticDiffEqAdaptiveAlgorithm ? 0.0 : 0.1
  srand(100)
  sol = solve(prob,alg;dt=dt)

  # Vector of element-wise absolute tolerances
  srand(100)
  sol2 = solve(prob,alg;dt=dt,abstol=fill(1e-2,4,2))

  @test sol.t == sol2.t && sol.u == sol2.u

  # Vector of element-wise relative tolerances
  srand(100)
  sol2 = solve(prob,alg;dt=dt,reltol=fill(1e-2,4,2))

  @test sol.t == sol2.t && sol.u == sol2.u
end
