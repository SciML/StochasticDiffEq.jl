using StochasticDiffEq, DiffEqDevTools, DiffEqProblemLibrary, Base.Test
srand(100)
prob = prob_sde_linear

## Solve and plot
println("Solve and Plot")
sol = solve(prob,EM(),dt=1//2^(4),tstops = [0.33])

@test 0.33 ∈ sol.t

sol = solve(prob,SRIW1(),tstops = [0.33])

@test 0.33 ∈ sol.t
