using StochasticDiffEq, DiffEqDevTools, DiffEqBase, Test
using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems; importsdeproblems()
import DiffEqProblemLibrary.SDEProblemLibrary: prob_sde_linear
Random.seed!(100)
prob = prob_sde_linear

integrator = init(prob,EM(),dt=1//2^(4),tstops = [0.33])

for (u,t) in tuples(integrator)
  @show u,t
end

sol = solve(prob,EM(),dt=1//2^(4),tstops = [0.33])

@test 0.33 ∈ sol.t

sol = solve(prob,EM(),tstops = [0.33,0.80,1.0])

@test sol.t == [0.0,0.33,0.80,1.0]

sol = solve(prob,SRIW1(),tstops = [0.33])

@test 0.33 ∈ sol.t
