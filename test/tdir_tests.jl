using StochasticDiffEq, DiffEqDevTools, DiffEqBase, DiffEqProblemLibrary, Base.Test
srand(100)
prob = deepcopy(prob_sde_linear)
prob2 = SDEProblem(prob.f,prob.g,prob.u0,(1.0,0.0))

sol = solve(prob2,EM(),dt=-1//2^(4),tstops = [0.33])

@test sol.t[15] > sol.t[16]

sol = solve(prob2,SRIW1(),tstops = [0.33])

@test sol.t[15] > sol.t[16]

prob = deepcopy(prob_sde_2Dlinear)
prob2 = SDEProblem(prob.f,prob.g,prob.u0,(1.0,0.0))

sol = solve(prob2,EM(),dt=-1//2^(4),tstops = [0.33])

@test sol.t[15] > sol.t[16]

sol = solve(prob2,SRIW1(),tstops = [0.33])

@test sol.t[15] > sol.t[16]
