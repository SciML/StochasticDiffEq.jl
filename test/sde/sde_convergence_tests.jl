using StochasticDiffEq, DiffEqDevTools, Test, Random
using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems; importsdeproblems()
using DiffEqProblemLibrary.SDEProblemLibrary: prob_sde_wave, prob_sde_cubic, prob_sde_additive
Random.seed!(100)
dts = (1/2) .^ (10:-1:2) #14->7 good plot

print("prob_sde_wave")
prob = prob_sde_wave
sim  = test_convergence(dts,prob,ImplicitEM(),numMonte=Int(1e1))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,ImplicitEM(nlsolve=StochasticDiffEq.NLFunctional()),numMonte=Int(1e1))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,ImplicitRKMil(),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim  = test_convergence(dts,prob,EM(),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,ISSEM(),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,ISSEM(nlsolve=StochasticDiffEq.NLFunctional()),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,LambaEM(),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim2 = test_convergence(dts,prob,RKMil(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
print(".")
sim2 = test_convergence(dts,prob,WangLi3SMil_A(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_B(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_C(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_D(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_E(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_F(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
print(".")
sim2 = test_convergence(dts,prob,SROCK1(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,SROCK2(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,SROCKEM(strong_order_1=false),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-0.5) < 0.2
sim2 = test_convergence(dts,prob,SROCKEM(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,SKSROCK(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-0.5) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=1),numMonte=Int(2e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=2),numMonte=Int(2e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=3),numMonte=Int(2e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=4),numMonte=Int(2e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=5),numMonte=Int(2e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim3 = test_convergence(dts,prob,SRI(),numMonte=Int(1e1))
@test abs(sim3.𝒪est[:final]-1.5) < 0.3
sim4 = test_convergence(dts,prob,SRIW1(),numMonte=Int(1e1))
@test abs(sim4.𝒪est[:final]-1.5) < 0.3
sim5 = test_convergence(dts,prob,SRIW2(),numMonte=Int(1e1))
@test abs(sim5.𝒪est[:final]-1.5) < 0.3
sim6 = test_convergence(dts,prob,SOSRI(),numMonte=Int(1e1))
@test abs(sim6.𝒪est[:final]-1.5) < 0.3
sim7 = test_convergence(dts,prob,SOSRI2(),numMonte=Int(1e1))
@test abs(sim7.𝒪est[:final]-1.5) < 0.3
println()

print("prob_sde_cubic")
prob = prob_sde_cubic
sim  = test_convergence(dts,prob,EM(),numMonte=Int(1e1))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,LambaEM(),numMonte=Int(1e1))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,ISSEM(),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,ImplicitEM(),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l2]-.5) < 0.2
sim  = test_convergence(dts,prob,ImplicitRKMil(),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim2 = test_convergence(dts,prob,RKMil(),numMonte=Int(1e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.22
print(".")
sim2 = test_convergence(dts,prob,WangLi3SMil_A(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_B(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_C(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_D(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_E(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_F(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
print(".")
sim2 = test_convergence(dts,prob,SROCK1(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,SROCK2(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,SROCKEM(strong_order_1=false),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-0.5) < 0.2
sim2 = test_convergence(dts,prob,SROCKEM(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,SKSROCK(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-0.5) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=1),numMonte=Int(2e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=2),numMonte=Int(2e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=3),numMonte=Int(2e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=4),numMonte=Int(2e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=5),numMonte=Int(2e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim3 = test_convergence(dts,prob,SRI(),numMonte=Int(1e1))
@test abs(sim3.𝒪est[:final]-1.5) < 0.3
sim4 = test_convergence(dts,prob,SRIW1(),numMonte=Int(1e1))
@test abs(sim4.𝒪est[:final]-1.5) < 0.3
sim5 = test_convergence(dts,prob,SRIW2(),numMonte=Int(1e1))
@test abs(sim5.𝒪est[:final]-1.5) < 0.3
sim6 = test_convergence(dts,prob,SOSRI(),numMonte=Int(1e1))
@test abs(sim6.𝒪est[:final]-1.5) < 0.3
sim7 = test_convergence(dts,prob,SOSRI2(),numMonte=Int(1e1))
@test abs(sim7.𝒪est[:final]-1.5) < 0.3
println()

print("prob_sde_additive")
prob = prob_sde_additive
sim  = test_convergence(dts,prob,EM(),numMonte=Int(1e1))
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim  = test_convergence(dts,prob,LambaEM(),numMonte=Int(1e1))
@test abs(sim.𝒪est[:l2]-1) < 0.2
dts = (1/2) .^ (10:-1:1)
sim  = test_convergence(dts,prob,ISSEM(),numMonte=Int(1e3))
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim  = test_convergence(dts,prob,ImplicitEM(),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim  = test_convergence(dts,prob,ImplicitRKMil(),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim2 = test_convergence(dts,prob,RKMil(),numMonte=Int(1e1))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
print(".")
sim2 = test_convergence(dts,prob,WangLi3SMil_A(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_B(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_C(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_D(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_E(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,WangLi3SMil_F(),numMonte=Int(2e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
print(".")
sim2 = test_convergence(dts,prob,SROCK1(),numMonte=Int(1e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,SROCK2(),numMonte=Int(1e2))
@test abs(sim2.𝒪est[:l∞]-2) < 0.2
@time sim2 = test_convergence(dts,prob,SROCKEM(strong_order_1=false),numMonte=Int(1e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
@time sim2 = test_convergence(dts,prob,SROCKEM(),numMonte=Int(1e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim2 = test_convergence(dts,prob,SKSROCK(),numMonte=Int(1e2))
@test abs(sim2.𝒪est[:l∞]-1) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=1),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=2),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=3),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=4),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=5),numMonte=Int(1e2))
@test abs(sim.𝒪est[:l∞]- 1) < 0.2
sim3 = test_convergence(dts,prob,SRI(),numMonte=Int(1e1))
@test abs(sim3.𝒪est[:final]-2) < 0.3
sim3 = test_convergence(dts,prob,SRIW2(),numMonte=Int(1e1))
@test abs(sim3.𝒪est[:final]-2) < 0.3

sim3 = test_convergence(dts,prob,SOSRI(),numMonte=Int(1e1))
@test abs(sim3.𝒪est[:final]-2) < 0.3
sim3 = test_convergence(dts,prob,SOSRI2(),numMonte=Int(1e1))
@test abs(sim3.𝒪est[:final]-2) < 0.3
print(".")
sim4 = test_convergence(dts,prob,SRA(),numMonte=Int(1e1))
@test abs(sim4.𝒪est[:final]-2) < 0.3
sim5 = test_convergence(dts,prob,SRA1(),numMonte=Int(1e1))
@test abs(sim5.𝒪est[:final]-2) < 0.3
sim6 = test_convergence(dts,prob,SRA2(),numMonte=Int(1e1))
@test abs(sim6.𝒪est[:final]-2) < 0.3
sim7 = test_convergence(dts,prob,SRA3(),numMonte=Int(1e1))
@test abs(sim7.𝒪est[:final]-2.5) < 0.3
sim8 = test_convergence(dts,prob,SOSRA(),numMonte=Int(1e1))
@test abs(sim8.𝒪est[:final]-2) < 0.3
sim9 = test_convergence(dts,prob,SOSRA2(),numMonte=Int(1e1))
@test abs(sim9.𝒪est[:final]-2) < 0.3
print(".")
dts = (1/2) .^ (10:-1:5) #14->7 good plot
sim10 = test_convergence(dts,prob,SKenCarp(),numMonte=Int(1e1))
@test abs(sim10.𝒪est[:final]-2) < 0.3
sim10 = test_convergence(dts,prob,SKenCarp(nlsolve=StochasticDiffEq.NLFunctional()),numMonte=Int(1e1))
@test abs(sim10.𝒪est[:final]-2) < 0.3

sim2 = test_convergence(dts,prob,SRA(tableau=StochasticDiffEq.constructSRA2()),numMonte=Int(1e1))
@test abs(sim2.𝒪est[:final]-2) < 0.3
sim3 = test_convergence(dts,prob,SRA(tableau=StochasticDiffEq.constructSRA3()),numMonte=Int(1e2))
@test abs(sim3.𝒪est[:final]-2.0) < 0.3
sim6 = test_convergence(dts,prob,SRIW1(),numMonte=Int(1e1))
@test abs(sim6.𝒪est[:final]-2) < 0.3
sim2 = test_convergence(dts,prob,SRA(tableau=StochasticDiffEq.constructExplicitSKenCarp()),numMonte=Int(1e1))
@test abs(sim2.𝒪est[:final]-2) < 0.3
println(".")
