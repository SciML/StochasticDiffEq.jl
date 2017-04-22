using StochasticDiffEq, DiffEqDevTools, DiffEqBase, DiffEqProblemLibrary, Base.Test
srand(100)
prob = deepcopy(prob_sde_linear)
prob.tspan = (1.0,0.0)

sol = solve(prob,EM(),dt=-1//2^(4),tstops = [0.33])

@test sol.t[15] > sol.t[16]

sol = solve(prob,SRIW1(),tstops = [0.33])

@test sol.t[15] > sol.t[16]

prob = deepcopy(prob_sde_2Dlinear)
prob.tspan = (1.0,0.0)

sol = solve(prob,EM(),dt=-1//2^(4),tstops = [0.33])

@test sol.t[15] > sol.t[16]

sol = solve(prob,SRIW1(),tstops = [0.33])

@test sol.t[15] > sol.t[16]
