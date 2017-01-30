using StochasticDiffEq, DiffEqDevTools, DiffEqProblemLibrary, Base.Test
prob = prob_sde_linear

srand(100)
sol1 = solve(prob,SRI(),dt=1//2^(4))
srand(100)
sol2 = solve(prob,SRIW1(),dt=1//2^(4))

@test sol1[end] ≈ sol2[end]

prob = prob_sde_2Dlinear

srand(100)
sol1 = solve(prob,SRI(),dt=1//2^(4))
srand(100)
sol2 = solve(prob,SRIW1(),dt=1//2^(4))

@test sol1[end] ≈ sol2[end]

prob = prob_sde_linear

srand(100)
integrator1 = init(prob,SRI(),dt=1//2^(4))
step!(integrator1); step!(integrator1)

srand(100)
integrator2 = init(prob,SRIW1(),dt=1//2^(4))
step!(integrator2); step!(integrator2)

@test integrator1.EEst ≈ integrator2.EEst

prob = prob_sde_2Dlinear

srand(100)
integrator1 = init(prob,SRI(),dt=1//2^(4))
step!(integrator1); step!(integrator1)

srand(100)
integrator2 = init(prob,SRIW1(),dt=1//2^(4))
step!(integrator2); step!(integrator2)

@test integrator1.EEst ≈ integrator2.EEst
