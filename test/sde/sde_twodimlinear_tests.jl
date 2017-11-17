using StochasticDiffEq, DiffEqProblemLibrary, DiffEqDevTools, Base.Test
srand(100)
prob = prob_sde_2Dlinear

## Solve and plot
println("Solve and Plot")
sol = solve(prob,EM(),dt=1/2^(3))
sol = solve(prob,RKMil(),dt=1/2^(3))
sol = solve(prob,SRI(),dt=1/2^(3))
sol = solve(prob,SRIW1(),dt=1/2^(3))
sol = solve(prob,ImplicitEM(),dt=1/2^(3))
sol = solve(prob,ImplicitEM(autodiff=false),dt=1/2^(3))
sol = solve(prob,ImplicitRKMil(),dt=1/2^(3))

sol = solve(prob,SRIW1(),dt=1/2^(3),save_everystep=false)


sol = solve(prob,SRIW1(),dt=1/2^(3),progress=true,progress_steps=1)

## Convergence Testing
println("Convergence Test on 2D Linear")
dts = 1./2.^(7:-1:4) #14->7 good plot

sim = test_convergence(dts,prob,EM(),numMonte=1000)
@test abs(sim.𝒪est[:l2]-.5) < 0.1

sim = test_convergence(dts,prob,ImplicitEM(),numMonte=100)
@test abs(sim.𝒪est[:l2]-.5) < 0.1

sim = test_convergence(dts,prob,ImplicitEM(theta=1),numMonte=100)
@test abs(sim.𝒪est[:l2]-.5) < 0.1

sim = test_convergence(dts,prob,ImplicitEM(symplectic=true),numMonte=100)
@test abs(sim.𝒪est[:l2]-.5) < 0.1

sim = test_convergence(dts,prob,ImplicitEM(symplectic=true,autodiff=false),numMonte=100)
@test abs(sim.𝒪est[:l2]-.5) < 0.1

sim = test_convergence(dts,prob,ImplicitRKMil(),numMonte=100)
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim = test_convergence(dts,prob,ImplicitRKMil(theta=1),numMonte=100)
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim = test_convergence(dts,prob,ImplicitRKMil(theta=1,autodiff=false),numMonte=200)
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim = test_convergence(dts,prob,ImplicitRKMil(symplectic=true),numMonte=100)
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim = test_convergence(dts,prob,ImplicitRKMil(symplectic=true,autodiff=false),numMonte=100)
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim2 = test_convergence(dts,prob,RKMil(),numMonte=100)
@test abs(sim2.𝒪est[:l∞]-1) < 0.2

sim3 = test_convergence(dts,prob,SRI(),numMonte=10)
@test abs(sim3.𝒪est[:final]-1.5) < 0.3

sim4 = test_convergence(dts,prob,SRIW1(),numMonte=100)
@test abs(sim4.𝒪est[:final]-1.5) < 0.3
