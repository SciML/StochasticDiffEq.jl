using StochasticDiffEq, Test, Random
using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems; importsdeproblems()
import DiffEqProblemLibrary.SDEProblemLibrary: prob_sde_stiffquadito

Random.seed!(100)
prob = prob_sde_stiffquadito
prob = remake(prob;p=(1e5,2.))
alg = AutoSOSRA2(SKenCarp(), maxstiffstep=5, maxnonstiffstep=2, stiffalgfirst=true)
@test StochasticDiffEq.isadaptive(alg)
@time sol = solve(prob, alg);
@test typeof(alg.algs[sol.alg_choice[1]]) <: SKenCarp
@test length(unique(sol.alg_choice)) == 2

Random.seed!(100)
prob = prob_sde_stiffquadito
prob = remake(prob;p=(1e5,2.))
@time sol = solve(prob, AutoSOSRI2(ImplicitRKMil(), maxstiffstep=2, maxnonstiffstep=2, stiffalgfirst=true),
                        abstol=1e-4,reltol=1e-4);
@test length(unique(sol.alg_choice)) == 2
