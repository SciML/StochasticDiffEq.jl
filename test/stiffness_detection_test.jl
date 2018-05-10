using StochasticDiffEq, DiffEqProblemLibrary, Base.Test

srand(100)
prob = DiffEqProblemLibrary.generate_stiff_quad(1e5,2.); # Not actually additive!
alg = AutoSOSRA2(SKenCarp(), maxstiffstep=5, maxnonstiffstep=2, stiffalgfirst=true)
@test StochasticDiffEq.isadaptive(alg)
@time sol = solve(prob, alg);
@test typeof(alg.algs[sol.alg_choice[1]]) <: SKenCarp
@test length(unique(sol.alg_choice)) == 2

srand(100)
prob = DiffEqProblemLibrary.generate_stiff_quad(1e5,2.);
@time sol = solve(prob, AutoSOSRI2(ImplicitRKMil(), maxstiffstep=3, maxnonstiffstep=3),
                        abstol=1e-4,reltol=1e-4);
@test length(unique(sol.alg_choice)) == 2
