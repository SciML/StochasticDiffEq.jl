using StochasticDiffEq, DiffEqProblemLibrary, Base.Test
prob = DiffEqProblemLibrary.generate_stiff_quad(1000.,2.)

srand(100)
@test StochasticDiffEq.isadaptive(AutoSOSRA2(SKenCarp()))
@time sol = solve(prob, AutoSOSRA2(SKenCarp()))
