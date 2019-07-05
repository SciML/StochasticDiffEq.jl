using StochasticDiffEq, Test, Random
using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems; importsdeproblems()
using DiffEqProblemLibrary.SDEProblemLibrary: prob_sde_linear, prob_sde_2Dlinear

Random.seed!(100)

for prob in (prob_sde_linear, prob_sde_2Dlinear)
  prob2 = remake(prob; tspan = (1.0, 0.0))

  solEM = solve(prob2, EM(), dt = -1//2^(4), tstops = [0.33])
  @test solEM.t[15] > solEM.t[16]

  solSRIW1 = solve(prob2, SRIW1(), tstops = [0.33])
  @test solSRIW1.t[15] > solSRIW1.t[16]
end
