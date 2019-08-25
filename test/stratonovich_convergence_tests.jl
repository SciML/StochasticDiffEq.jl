using StochasticDiffEq, Test, Random, DiffEqDevTools
using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems; importsdeproblems()
using DiffEqProblemLibrary.SDEProblemLibrary: prob_sde_linear_stratonovich, prob_sde_2Dlinear_stratonovich
Random.seed!(100)
dts = 1 ./2 .^(10:-1:2) #14->7 good plot

prob = prob_sde_linear_stratonovich
sim  = test_convergence(dts,prob,EulerHeun(),trajectories=Int(5e2))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,LambaEulerHeun(),trajectories=Int(5e2))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,ISSEulerHeun(),trajectories=Int(5e2))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,ImplicitEulerHeun(),trajectories=Int(5e2))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,ImplicitEulerHeun(theta=1),trajectories=Int(5e2))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,ImplicitEulerHeun(symplectic=true),trajectories=Int(5e2))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim = test_convergence(dts,prob,RKMil(interpretation=:Stratonovich),trajectories=Int(5e2))
@test abs(sim.𝒪est[:l2]-1) < 0.2

sim = test_convergence(dts,prob,RKMil_General(interpretation=:Stratonovich),trajectories=Int(5e2))
@test abs(sim.𝒪est[:l2]-1) < 0.2

sim  = test_convergence(dts,prob,ImplicitRKMil(interpretation=:Stratonovich),trajectories=Int(5e2))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,SROCK1(interpretation=:Stratonovich),trajectories=Int(2e2))
@test abs(sim.𝒪est[:l2]-1) < 0.15

sim  = test_convergence(dts,prob,KomBurSROCK2(),trajectories=Int(2e2))
@test abs(sim.𝒪est[:final]-2) < 0.20

println("Now 2D")

prob = prob_sde_2Dlinear_stratonovich

sim  = test_convergence(dts,prob,EulerHeun(),trajectories=Int(5e1))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,LambaEulerHeun(),trajectories=Int(5e1))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,ISSEulerHeun(),trajectories=Int(5e1))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,ImplicitEulerHeun(),trajectories=Int(5e1))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,ImplicitEulerHeun(theta=1),trajectories=Int(5e1))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,ImplicitEulerHeun(symplectic=true),trajectories=Int(5e1))
@test abs(sim.𝒪est[:l2]-1) < 0.1

println("RKMils")

sim = test_convergence(dts,prob,RKMil(interpretation=:Stratonovich),trajectories=Int(1e2))
@test abs(sim.𝒪est[:l2]-1) < 0.2

sim = test_convergence(dts,prob,RKMil_General(interpretation=:Stratonovich),trajectories=Int(1e2))
@test abs(sim.𝒪est[:l2]-1) < 0.2

sim  = test_convergence(dts,prob,ImplicitRKMil(interpretation=:Stratonovich),
                        trajectories=Int(1e2))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,ImplicitRKMil(theta=1,interpretation=:Stratonovich),
                        trajectories=Int(1e2))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,ImplicitRKMil(symplectic=true,interpretation=:Stratonovich),
                        trajectories=Int(1e2))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,SROCK1(interpretation=:Stratonovich),trajectories=Int(1e2))
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,KomBurSROCK2(),trajectories=Int(2e2))
@test abs(sim.𝒪est[:final]-2) < 0.20

Random.seed!(200)
sol = solve(prob,EulerHeun(),dt=1/4)
