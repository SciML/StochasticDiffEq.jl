using StochasticDiffEq, DiffEqDevTools, Test, Random
using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems; importsdeproblems()
using DiffEqProblemLibrary.SDEProblemLibrary: prob_sde_wave, prob_sde_cubic
Random.seed!(100)
dts = (1/2) .^ (10:-1:2) #14->7 good plot

print("prob_sde_wave")
prob = prob_sde_wave
sim  = test_convergence(dts,prob,ImplicitEM(),trajectories=Int(1e1))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,ImplicitEM(nlsolve=StochasticDiffEq.NLFunctional()),trajectories=Int(1e1))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,ImplicitRKMil(),trajectories=Int(1e2))
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim  = test_convergence(dts,prob,EM(),trajectories=Int(1e2))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,ISSEM(),trajectories=Int(1e2))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,ISSEM(nlsolve=StochasticDiffEq.NLFunctional()),trajectories=Int(1e2))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,LambaEM(),trajectories=Int(1e2))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim2 = test_convergence(dts,prob,RKMil(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
print(".")
sim2 = test_convergence(dts,prob,WangLi3SMil_A(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_B(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_C(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_D(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_E(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_F(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
print(".")
sim2 = test_convergence(dts,prob,SROCK1(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,SROCK2(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,SROCKEM(strong_order_1=false),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-0.5) < 0.2
sim2 = test_convergence(dts,prob,SROCKEM(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,SKSROCK(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-0.5) < 0.2
sim2 = test_convergence(dts,prob,SROCKC2(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2

#omitting tests for incomplete methods
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=1),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=2),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=3),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=4),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=5),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2

sim3 = test_convergence(dts,prob,SRI(),trajectories=Int(1e1))
@test abs(sim3.𝒪est[:final]-1.5) < 0.3
sim4 = test_convergence(dts,prob,SRIW1(),trajectories=Int(1e1))
@test abs(sim4.𝒪est[:final]-1.5) < 0.3
sim5 = test_convergence(dts,prob,SRIW2(),trajectories=Int(1e1))
@test abs(sim5.𝒪est[:final]-1.5) < 0.3
sim6 = test_convergence(dts,prob,SOSRI(),trajectories=Int(1e1))
@test abs(sim6.𝒪est[:final]-1.5) < 0.3
sim7 = test_convergence(dts,prob,SOSRI2(),trajectories=Int(1e1))
@test abs(sim7.𝒪est[:final]-1.5) < 0.3
println()

print("prob_sde_cubic")
prob = prob_sde_cubic
sim  = test_convergence(dts,prob,EM(),trajectories=Int(1e1))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,LambaEM(),trajectories=Int(1e1))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,ISSEM(),trajectories=Int(1e2))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,ImplicitEM(),trajectories=Int(1e2))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,ImplicitRKMil(),trajectories=Int(1e2))
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim2 = test_convergence(dts,prob,RKMil(),trajectories=Int(1e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.22
print(".")
sim2 = test_convergence(dts,prob,WangLi3SMil_A(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_B(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_C(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_D(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_E(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_F(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
print(".")
sim2 = test_convergence(dts,prob,SROCK1(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,SROCK2(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,SROCKEM(strong_order_1=false),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-0.5) < 0.2
sim2 = test_convergence(dts,prob,SROCKEM(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,SKSROCK(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-0.5) < 0.2
sim2 = test_convergence(dts,prob,SROCKC2(),trajectories=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2

#omitting tests for incomplete methods
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=1),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=2),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=3),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=4),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=5),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2

sim3 = test_convergence(dts,prob,SRI(),trajectories=Int(1e1))
@test abs(sim3.𝒪est[:final]-1.5) < 0.3
sim4 = test_convergence(dts,prob,SRIW1(),trajectories=Int(1e1))
@test abs(sim4.𝒪est[:final]-1.5) < 0.3
sim5 = test_convergence(dts,prob,SRIW2(),trajectories=Int(1e1))
@test abs(sim5.𝒪est[:final]-1.5) < 0.3
sim6 = test_convergence(dts,prob,SOSRI(),trajectories=Int(1e1))
@test abs(sim6.𝒪est[:final]-1.5) < 0.3
sim7 = test_convergence(dts,prob,SOSRI2(),trajectories=Int(1e1))
@test abs(sim7.𝒪est[:final]-1.5) < 0.3
println()
