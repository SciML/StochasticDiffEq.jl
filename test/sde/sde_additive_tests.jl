using StochasticDiffEq
srand(100)

prob = prob_sde_additive
sol =solve(prob,SRA(),dt=1/2^(3))
sol =solve(prob,SRA1(),dt=1/2^(3))

prob = prob_sde_additivesystem

## Solve and plot
println("Solve and Plot")
sol =solve(prob,SRA(),dt=1/2^(3))
sol =solve(prob,SRA1(),dt=1/2^(3))

#Now do the simulation 10000 times in parallel. Return an array
solArr = monte_carlo_simulation(prob,SRIW1(),dt=1//2^(3),numMonte=5)

#First index is the sime, so sol.timeseries[1,..] is the initial condition
#Last indices are the indexes of the variables. Since our initial condition
#Has 4 rows and two columns, sol.timeseries[..,1] returns the time series for the
#first row, and sol.timeseries[..,2] returns the time series for the second.
TEST_PLOT && plot(sol,plot_analytic=true)

## Convergence Testing
println("Convergence Test on MultiDimAdditive")
dts = 1./2.^(7:-1:4) #14->7 good plot

sim = test_convergence(dts,prob,SRA(),numMonte=5)

sim2 = test_convergence(dts,prob,SRA1(),numMonte=5)

@test abs(sim.𝒪est[:l2]-2) + abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low dts for testing!
