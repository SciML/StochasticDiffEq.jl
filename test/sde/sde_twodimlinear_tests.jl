using StochasticDiffEq, DiffEqProblemLibrary, DiffEqDevTools
srand(70)
prob = prob_sde_2Dlinear

## Solve and plot
println("Solve and Plot")
sol = solve(prob,EM(),dt=1/2^(3))
sol = solve(prob,RKMil(),dt=1/2^(3))
sol = solve(prob,SRI(),dt=1/2^(3))
sol = solve(prob,SRIW1(),dt=1/2^(3))

sol = solve(prob,SRIW1(),dt=1/2^(3),save_timeseries=false)


sol = solve(prob,SRIW1(),dt=1/2^(3),progressbar=true,progress_steps=1)


#Now do the simulation 5 times in parallel. Return an array
solArr = monte_carlo_simulation(prob,SRIW1(),dt=1//2^(3),numMonte=5)

TEST_PLOT && plot(sol,plot_analytic=true)

## Convergence Testing
println("Convergence Test on 2D Linear")
dts = 1./2.^(7:-1:4) #14->7 good plot

sim = test_convergence(dts,prob,EM(),numMonte=10)
@test abs(sim.𝒪est[:l2]-.5) < 0.1

sim2 = test_convergence(dts,prob,RKMil(),numMonte=10)
@test abs(sim2.𝒪est[:l∞]-1) < 0.1

sim3 = test_convergence(dts,prob,SRI(),numMonte=10)
@test abs(sim3.𝒪est[:final]-1.5) < 0.3

sim4 = test_convergence(dts,prob,SRIW1(),numMonte=10,save_timeseries=false)
@test abs(sim4.𝒪est[:final]-1.5) < 0.3
