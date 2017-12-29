using StochasticDiffEq, DiffEqProblemLibrary, DiffEqDevTools, Base.Test
srand(100)

f(t,u) = 0.0
f(::Type{Val{:analytic}},t,u0,W) = W
g(t,u) = 1.0
prob = SDEProblem(f,g,0.0,(0.0,1.0))
sol1 = solve(prob,SRA1(),dt=1/2^(3),seed=1,adaptive=false)
sol2 = solve(prob,SOSRA(),dt=1/2^(3),seed=1,adaptive=false)
sol3 = solve(prob,RackKenCarp(),dt=1/2^(3),seed=1,adaptive=false)
@test sol1.errors[:l2] ≈ 0.0 atol = 1e-14
@test sol2.errors[:l2] ≈ 0.0 atol = 1e-14
@test sol3.errors[:l2] ≈ 0.0 atol = 1e-14

prob = prob_sde_additive
sol =solve(prob,SRA(),dt=1/2^(3))
sol =solve(prob,SRA1(),dt=1/2^(3))
sol =solve(prob,SRA2(),dt=1/2^(3))
sol =solve(prob,SRA3(),dt=1/2^(3))
sol =solve(prob,SOSRA(),dt=1/2^(3))
sol =solve(prob,SOSRA2(),dt=1/2^(3))
sol =solve(prob,RackKenCarp(),dt=1/2^(3))

prob = prob_sde_additivesystem

sol =solve(prob,SRA(),dt=1/2^(3))
sol =solve(prob,SRA1(),dt=1/2^(3))
sol =solve(prob,SRA2(),dt=1/2^(3))
sol =solve(prob,SRA3(),dt=1/2^(3))
sol =solve(prob,SOSRA(),dt=1/2^(3))
sol =solve(prob,SOSRA2(),dt=1/2^(3))
sol =solve(prob,RackKenCarp(),dt=1/2^(3))

## Convergence Testing
println("Convergence Test on MultiDimAdditive")
dts = 1./2.^(7:-1:4) #14->7 good plot

sim = test_convergence(dts,prob,SRA(),numMonte=5)
@test abs(sim.𝒪est[:l2]-2) < 0.1
sim2 = test_convergence(dts,prob,SRA1(),numMonte=5)
@test abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low dts for testing!
sim2 = test_convergence(dts,prob,SRA2(),numMonte=5)
@test abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low dts for testing!
sim2 = test_convergence(dts,prob,SRA3(),numMonte=5)
@test abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low dts for testing!
sim2 = test_convergence(dts,prob,SOSRA(),numMonte=10)
@test abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low dts for testing!
sim2 = test_convergence(dts,prob,SOSRA2(),numMonte=5)
@test abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low dts for testing!
sim2 = test_convergence(dts,prob,RackKenCarp(),numMonte=20)
@test abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low dts for testing!
