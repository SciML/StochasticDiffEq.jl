using StochasticDiffEq, DiffEqDevTools, Test, Random
using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems; importsdeproblems()
import DiffEqProblemLibrary.SDEProblemLibrary: prob_sde_additive, prob_sde_additivesystem
Random.seed!(100)

println("Bunch of additive solves")

f_bm(u,p,t) = 0.0
f_bm(::Type{Val{:analytic}},u0,p,t,W) = W
g_bm(u,p,t) = 1.0
prob = SDEProblem(f_bm,g_bm,0.0,(0.0,1.0))
sol1 = solve(prob,SRA1(),dt=1/2^(3),adaptive=false)
sol2 = solve(prob,SOSRA(),dt=1/2^(3),adaptive=false)
sol3 = solve(prob,SKenCarp(),dt=1/2^(3),adaptive=false)
@test sol1.errors[:l2] ≈ 0.0 atol = 1e-14
@test sol2.errors[:l2] ≈ 0.0 atol = 1e-14
@test sol3.errors[:l2] ≈ 0.0 atol = 1e-14

prob = prob_sde_additive

# Test error in stepping and seeding simultaniously
sol  = solve(prob,SRA(StochasticDiffEq.constructSOSRA()),dt=1/2^(3),seed=1)
sol2 = solve(prob,SOSRA(),dt=1/2^(3),seed=1)
@test sol.t ≈ sol2.t

sol =solve(prob,SRA1(),dt=1/2^(3))
sol =solve(prob,SRA2(),dt=1/2^(3))
sol =solve(prob,SRA3(),dt=1/2^(3))
sol =solve(prob,SOSRA(),dt=1/2^(3))
sol =solve(prob,SOSRA2(),dt=1/2^(3))
sol =solve(prob,SKenCarp(),dt=1/2^(3))
sol =solve(prob,SKenCarp(nlsolve=StochasticDiffEq.NLNewton()),dt=1/2^(3))

prob = prob_sde_additivesystem

sol =solve(prob,SRA(),dt=1/2^(3))
sol =solve(prob,SRA1(),dt=1/2^(3))
sol =solve(prob,SRA2(),dt=1/2^(3))
sol =solve(prob,SRA3(),dt=1/2^(3))
sol =solve(prob,SOSRA(),dt=1/2^(3))
sol =solve(prob,SOSRA2(),dt=1/2^(3))
sol =solve(prob,SKenCarp(),dt=1/2^(3))

## Convergence Testing
println("Convergence Test on MultiDimAdditive")
dts = (1/2) .^ (7:-1:4) #14->7 good plot

sim = test_convergence(dts,prob,SRA(),numMonte=10)
@test abs(sim.𝒪est[:l2]-2) < 0.1
sim2 = test_convergence(dts,prob,SRA1(),numMonte=10)
@test abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low dts for testing!
sim2 = test_convergence(dts,prob,SRA2(),numMonte=10)
@test abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low dts for testing!
sim2 = test_convergence(dts,prob,SRA3(),numMonte=10)
@test abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low dts for testing!
sim2 = test_convergence(dts,prob,SOSRA(),numMonte=10)
@test abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low dts for testing!
sim2 = test_convergence(dts,prob,SOSRA2(),numMonte=10)
@test abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low dts for testing!
dts = (1/2) .^ (14:-1:11) #14->7 good plot
Random.seed!(100)
sim2 = test_convergence(dts,prob,SKenCarp(),numMonte=20)
@test abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low dts for testing!
sim2 = test_convergence(dts,prob,SKenCarp(nlsolve=StochasticDiffEq.NLFunctional()),numMonte=20)
@test abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low dts for testing!
