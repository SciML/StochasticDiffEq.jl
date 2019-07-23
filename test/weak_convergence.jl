using Distributed
@everywhere using StochasticDiffEq, DiffEqDevTools, Test
@everywhere using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems
@everywhere importsdeproblems()
using DiffEqProblemLibrary.SDEProblemLibrary: prob_sde_linear, prob_sde_2Dlinear, prob_sde_additive
Random.seed!(100)
dts = 1./2.^(10:-1:2) #14->7 good plot

prob = prob_sde_linear
sim  = test_convergence(dts,prob,EM(),trajectories=Int(1e4),
                        weak_timeseries_errors=true,weak_dense_errors=true)
@test abs(sim.𝒪est[:weak_final]-1) < 0.3
@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
@test abs(sim.𝒪est[:weak_L2]-1) < 0.3
@test abs(sim.𝒪est[:weak_L∞]-1) < 0.3
sim2 = test_convergence(dts,prob,RKMil(),trajectories=Int(1e4),
                        weak_timeseries_errors=true,dense_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SROCK1(),trajectories=Int(1e4),
                        weak_timeseries_errors=true,dense_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SROCK2(),trajectories=Int(1e5),
                        weak_timeseries_errors=true,dense_errors=true)
@test abs(sim2.𝒪est[:weak_final]-2) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
sim2 = test_convergence(dts,prob,SROCKEM(strong_order_1=false),trajectories=Int(1e4),
                        weak_timeseries_errors=true,dense_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SROCKEM(),trajectories=Int(1e4),
                        weak_timeseries_errors=true,dense_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SKSROCK(),trajectories=Int(1e4),
                        weak_timeseries_errors=true,dense_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SROCKC2(),trajectories=Int(1e5),
                        weak_timeseries_errors=true,dense_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3

#omitting tests for incomplete methods
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=1),trajectories=Int(1e5),
#                           weak_timeseries_errors=true,dense_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=2),trajectories=Int(1e5),
#                           weak_timeseries_errors=true,dense_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=3),trajectories=Int(1e5),
#                           weak_timeseries_errors=true,dense_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=4),trajectories=Int(1e5),
#                           weak_timeseries_errors=true,dense_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=5),trajectories=Int(1e5),
#                           weak_timeseries_errors=true,dense_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3

sim2 = test_convergence(dts,prob,WangLi3SMil_A(),trajectories=Int(1e4),
                        weak_timeseries_errors=true,dense_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_B(),trajectories=Int(1e4),
                        weak_timeseries_errors=true,dense_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_C(),trajectories=Int(1e4),
                        weak_timeseries_errors=true,dense_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_D(),trajectories=Int(1e4),
                        weak_timeseries_errors=true,dense_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_E(),trajectories=Int(1e4),
                        weak_timeseries_errors=true,dense_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_F(),trajectories=Int(1e4),
                        weak_timeseries_errors=true,dense_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim3 = test_convergence(dts,prob,SRI(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim3.𝒪est[:weak_final]-2) < 0.3
@test abs(sim3.𝒪est[:weak_l2]-2) < 0.3
@test abs(sim3.𝒪est[:weak_l∞]-2) < 0.3
sim4 = test_convergence(dts,prob,SRIW1(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim4.𝒪est[:weak_final]-2) < 0.3
@test abs(sim4.𝒪est[:weak_l2]-2) < 0.3
@test abs(sim4.𝒪est[:weak_l∞]-2) < 0.3

prob = prob_sde_2Dlinear
sim  = test_convergence(dts,prob,EM(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim.𝒪est[:weak_final]-1) < 0.3
@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,RKMil(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SROCK1(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SROCK2(),trajectories=Int(1e5),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-2) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
sim2 = test_convergence(dts,prob,SROCKEM(strong_order_1=false),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SROCKEM(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SKSROCK(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SROCKC2(),trajectories=Int(1e5),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3

#omitting tests for incomplete methods
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=1),trajectories=Int(1e5),
#                    weak_timeseries_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=2),trajectories=Int(1e5),
#                    weak_timeseries_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=3),trajectories=Int(1e5),
#                    weak_timeseries_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=4),trajectories=Int(1e5),
#                    weak_timeseries_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=5),trajectories=Int(1e5),
#                    weak_timeseries_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3

sim2 = test_convergence(dts,prob,WangLi3SMil_A(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_B(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_C(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_D(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_E(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_F(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim3 = test_convergence(dts,prob,SRI(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim3.𝒪est[:weak_final]-2) < 0.3
@test abs(sim3.𝒪est[:weak_l2]-2) < 0.3
@test abs(sim3.𝒪est[:weak_l∞]-2) < 0.3
sim4 = test_convergence(dts,prob,SRIW1(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim4.𝒪est[:weak_final]-2) < 0.3
@test abs(sim4.𝒪est[:weak_l2]-2) < 0.3
@test abs(sim4.𝒪est[:weak_l∞]-2) < 0.35

prob = prob_sde_additive
sim  = test_convergence(dts,prob,EM(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim.𝒪est[:weak_final]-1) < 0.3
@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,RKMil(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SROCK1(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SROCK2(),trajectories=Int(1e5),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-2) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
sim2 = test_convergence(dts,prob,SROCKEM(strong_order_1=false),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SROCKEM(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SKSROCK(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,SROCKC2(),trajectories=Int(1e5),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3

#omitting tests for incomplete methods
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=1),trajectories=Int(1e5),
#                    weak_timeseries_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=2),trajectories=Int(1e5),
#                    weak_timeseries_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=3),trajectories=Int(1e5),
#                    weak_timeseries_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=4),trajectories=Int(1e5),
#                    weak_timeseries_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3
# sim2 = test_convergence(dts,prob,TangXiaoSROCK2(version_num=5),trajectories=Int(1e5),
#                    weak_timeseries_errors=true)
# @test abs(sim2.𝒪est[:weak_final]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l2]-2) < 0.3
# @test abs(sim2.𝒪est[:weak_l∞]-2) < 0.3

sim2 = test_convergence(dts,prob,WangLi3SMil_A(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_B(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_C(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_D(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_E(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim2 = test_convergence(dts,prob,WangLi3SMil_F(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim2.𝒪est[:weak_final]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l2]-1) < 0.3
@test abs(sim2.𝒪est[:weak_l∞]-1) < 0.3
sim3 = test_convergence(dts,prob,SRI(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim3.𝒪est[:weak_final]-2) < 0.3
@test abs(sim3.𝒪est[:weak_l2]-2) < 0.3
@test abs(sim3.𝒪est[:weak_l∞]-2) < 0.3
sim4 = test_convergence(dts,prob,SRIW1(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim4.𝒪est[:weak_final]-2) < 0.3
@test abs(sim4.𝒪est[:weak_l2]-2) < 0.3
@test abs(sim4.𝒪est[:weak_l∞]-2) < 0.35
sim5 = test_convergence(dts,prob,SRA(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim5.𝒪est[:weak_final]-2) < 0.3
@test abs(sim5.𝒪est[:weak_l2]-2) < 0.3
@test abs(sim5.𝒪est[:weak_l∞]-2) < 0.3
sim6 = test_convergence(dts,prob,SRA1(),trajectories=Int(1e4),
                        weak_timeseries_errors=true)
@test abs(sim6.𝒪est[:weak_final]-2) < 0.3
@test abs(sim6.𝒪est[:weak_l2]-2) < 0.3
@test abs(sim6.𝒪est[:weak_l∞]-2) < 0.3
