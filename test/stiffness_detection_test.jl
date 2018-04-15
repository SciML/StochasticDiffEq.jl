using StochasticDiffEq, DiffEqProblemLibrary, Base.Test
prob = DiffEqProblemLibrary.generate_stiff_quad(1e8,2.);

srand(100)
alg = AutoSOSRA2(SKenCarp(), maxstiffstep=2, maxnonstiffstep=2, stiffalgfirst=true)
@test StochasticDiffEq.isadaptive(alg)
@time sol = solve(prob, alg);
@test typeof(alg.algs[sol.alg_choice[1]]) <: SKenCarp
@test length(unique(sol.alg_choice)) == 2

srand(100)
@time sol = solve(prob, AutoSOSRI2(SKenCarp(), maxstiffstep=3, maxnonstiffstep=3));
@test length(unique(sol.alg_choice)) == 2
