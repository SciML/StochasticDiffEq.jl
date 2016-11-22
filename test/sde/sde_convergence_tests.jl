@everywhere using StochasticDiffEq
srand(100)
dts = 1./2.^(10:-1:4) #14->7 good plot

prob = prob_sde_wave
sim  = test_convergence(dts,prob,EM(),numMonte=Int(1e1))
sim2 = test_convergence(dts,prob,RKMil(),numMonte=Int(1e1))
sim3 = test_convergence(dts,prob,SRI(),numMonte=Int(1e1))
sim4 = test_convergence(dts,prob,SRIW1(),numMonte=Int(1e1))
sim5 = test_convergence(dts,prob,SRIVectorized(),numMonte=Int(1e1))

@test abs(sim.𝒪est[:l2]-.5) + abs(sim2.𝒪est[:l∞]-1) + abs(sim3.𝒪est[:final]-1.5) + abs(sim4.𝒪est[:final]-1.5) + abs(sim5.𝒪est[:final]-1.5) <.5 #High tolerance since low dts for testing!

prob = prob_sde_cubic
sim  = test_convergence(dts,prob,EM(),numMonte=Int(1e1))
sim2 = test_convergence(dts,prob,RKMil(),numMonte=Int(1e1))
sim3 = test_convergence(dts,prob,SRI(),numMonte=Int(1e1))
sim4 = test_convergence(dts,prob,SRIW1(),numMonte=Int(1e1))
@test abs(sim.𝒪est[:l2]-.5) + abs(sim2.𝒪est[:l∞]-1) + abs(sim3.𝒪est[:final]-1.5) + abs(sim4.𝒪est[:final]-1.5) <.6 #High tolerance since low dts for testing!

## Convergence Testing
prob = prob_sde_additive
sim  = test_convergence(dts,prob,EM(),numMonte=Int(1e1))
sim2 = test_convergence(dts,prob,RKMil(),numMonte=Int(1e1))
sim3 = test_convergence(dts,prob,SRI(),numMonte=Int(1e1))
sim4 = test_convergence(dts,prob,SRA(),numMonte=Int(1e1))
sim5 = test_convergence(dts,prob,SRA1(),numMonte=Int(1e1))
sim6 = test_convergence(dts,prob,SRIW1(),numMonte=Int(1e1))
sim7 = test_convergence(dts,prob,SRAVectorized(),numMonte=Int(1e1))
@test abs(sim.𝒪est[:l2]-1) + abs(sim2.𝒪est[:l∞]-1) + abs(sim3.𝒪est[:final]-2) + abs(sim4.𝒪est[:final]-2) + abs(sim5.𝒪est[:final]-2) + abs(sim6.𝒪est[:final]-2) + abs(sim7.𝒪est[:final]-2)  <.4 #High tolerance since low dts for testing!
