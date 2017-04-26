@everywhere using StochasticDiffEq, DiffEqProblemLibrary, DiffEqDevTools, Base.Test
srand(100)
dts = 1./2.^(10:-1:2) #14->7 good plot

sim = test_convergence(dts,prob_sde_wave,SRI(tableau=StochasticDiffEq.constructSRIOpt1()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-1.5) < 0.3
sim = test_convergence(dts,prob_sde_cubic,SRI(tableau=StochasticDiffEq.constructSRIOpt1()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-1.5) < 0.3
sim = test_convergence(dts,prob_sde_additive,SRI(tableau=StochasticDiffEq.constructSRIOpt1()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-2.0) < 0.3

sim = test_convergence(dts,prob_sde_wave,SRI(tableau=StochasticDiffEq.constructSRIOpt2()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-1.5) < 0.3
sim = test_convergence(dts,prob_sde_cubic,SRI(tableau=StochasticDiffEq.constructSRIOpt2()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-1.5) < 0.3
sim = test_convergence(dts,prob_sde_additive,SRI(tableau=StochasticDiffEq.constructSRIOpt2()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-2.0) < 0.3

sim = test_convergence(dts,prob_sde_wave,SRI(tableau=StochasticDiffEq.constructSRIOpt3()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-1.5) < 0.3
sim = test_convergence(dts,prob_sde_cubic,SRI(tableau=StochasticDiffEq.constructSRIOpt3()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-1.5) < 0.3
sim = test_convergence(dts,prob_sde_additive,SRI(tableau=StochasticDiffEq.constructSRIOpt3()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-2.0) < 0.3

sim = test_convergence(dts,prob_sde_wave,SRI(tableau=StochasticDiffEq.constructSRIOpt4()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-1.5) < 0.3
sim = test_convergence(dts,prob_sde_cubic,SRI(tableau=StochasticDiffEq.constructSRIOpt4()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-1.5) < 0.3
sim = test_convergence(dts,prob_sde_additive,SRI(tableau=StochasticDiffEq.constructSRIOpt4()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-2.0) < 0.3

sim = test_convergence(dts,prob_sde_wave,SRI(tableau=StochasticDiffEq.constructSRIOpt7()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-1.5) < 0.3
sim = test_convergence(dts,prob_sde_cubic,SRI(tableau=StochasticDiffEq.constructSRIOpt7()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-1.5) < 0.3
sim = test_convergence(dts,prob_sde_additive,SRI(tableau=StochasticDiffEq.constructSRIOpt7()),numMonte=Int(1e3))
@test abs(sim.𝒪est[:final]-2.0) < 0.3
