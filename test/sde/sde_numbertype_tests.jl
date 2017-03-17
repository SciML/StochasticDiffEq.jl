using StochasticDiffEq
srand(100)
prob = prob_sde_2Dlinear

## Solve and plot
println("Solve and Plot")
sol =solve(prob,SRI(),dt=1//2^(4),abstol=1,reltol=0)
err1 = sol.errors[:final]


println("Solve and Plot")
sol =solve(prob,SRI(),dt=1//2^(4),abstol=1,reltol=0)
err1 = sol.errors[:final]

@test true
