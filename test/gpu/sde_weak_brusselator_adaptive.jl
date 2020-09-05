using DiffEqGPU, StochasticDiffEq, Test, DiffEqNoiseProcess
using CUDA
using Random

function brusselator_f!(du,u,p,t)
 @inbounds begin
     du[1] = (p[1]-1)*u[1]+p[1]*u[1]^2+(u[1]+1)^2*u[2]
     du[2] = -p[1]*u[1]-p[1]*u[1]^2-(u[1]+1)^2*u[2]
 end
 nothing
end

function scalar_noise!(du,u,p,t)
 @inbounds begin
     du[1] = p[2]*u[1]*(1+u[1])
     du[2] = -p[2]*u[1]*(1+u[1])
 end
 nothing
end

function prob_func(prob, i, repeat)
    Random.seed!(seeds[i])
    W = WienerProcess(0.0,0.0,0.0)
    remake(prob,noise=W)
end


# fix seeds
seed = 100
Random.seed!(seed)
numtraj= 100
seeds = rand(UInt, numtraj)
W = WienerProcess(0.0,0.0,0.0)

CUDA.allowscalar(false)
u0 = [-0.1f0,0.0f0]
tspan = (0.0f0,100.0f0)
p = [1.9f0,0.1f0]

prob = SDEProblem(brusselator_f!,scalar_noise!,u0,tspan,p, noise=W)
ensembleprob = EnsembleProblem(prob, prob_func = prob_func)

@info "Brusselator"

#Performance check with nvvp
# CUDAnative.CUDAdrv.@profile
# check either on CPU with EnsembleCPUArray() or on GPU with EnsembleGPUArray()
sol = @time solve(ensembleprob,DRI1(),EnsembleCPUArray(),trajectories=numtraj)
#sol = @time solve(ensembleprob,DRI1(),EnsembleGPUArray(),trajectories=numtraj)



# using Plots; plotly()  # or gr()
# using Plots; plot(sol,linealpha=0.6,color=:blue,vars=(0,1),title="Phase Space Plot")
# plot!(sol,linealpha=0.6,color=:red,vars=(0,2),title="Phase Space Plot")
# plot(sol,linealpha=0.6,color=:grey, vars=(1,2),title="Phase Space Plot")
#
# scatter(sol,linealpha=0.6,color=:blue,vars=(0,1),title="Phase Space Plot")
#
# summ = EnsembleSummary(sol,0.0f0:0.5f0:100f0)
# plot(summ,fillalpha=0.5)
#
#
# using DifferentialEquations.EnsembleAnalysis
# meansol = timeseries_steps_mean(sol)
# x1 = [x[1] for x in meansol.u]
# x2 = [x[2] for x in meansol.u]
# dts = []
# tmp1 = tspan[1]
# for tmp2 in meansol.t
#   global tmp1
#   push!(dts,tmp2-tmp1)
#   tmp1 = tmp2
# end
# #
# plot(x1,x2)
# plt = plot(dts)
# savefig(plt, "chosen_timesteps.png")



###
# oop not yet supported by DiffEqGPU
###
#
# function brusselator_f(u,p,t)
#  @inbounds begin
#    dx = (p[1]-1)*u[1]+p[1]*u[1]^2+(u[1]+1)^2*u[2]
#    dy = -p[1]*u[1]-p[1]*u[1]^2-(u[1]+1)^2*u[2]
#  end
#  return [dx,dy]
# end
#
# function scalar_noise(u,p,t)
#  @inbounds begin
#      dx = p[2]*u[1]*(1+u[1])
#      dy = -p[2]*u[1]*(1+u[1])
#  end
#  return [dx,dy]
# end
#
# function prob_func(prob, i, repeat)
#     Random.seed!(seeds[i])
#     W = WienerProcess(0.0,0.0,0.0)
#     remake(prob,noise=W)
# end
#
#
# # fix seeds
# seed = 100
# Random.seed!(seed)
# numtraj= 10
# seeds = rand(UInt, numtraj)
# W = WienerProcess(0.0,0.0,0.0)
#
# #CUDA.allowscalar(false)
# u0 = [-0.1f0,0.0f0]
# tspan = (0.0f0,100.0f0)
# p = [1.9f0,0.1f0]
#
# prob = SDEProblem(brusselator_f,scalar_noise,u0,tspan,p, noise=W)
#
#
# ensembleprob = EnsembleProblem(prob, prob_func = prob_func)
# sol = solve(ensembleprob,DRI1(),EnsembleCPUArray(),trajectories=numtraj)
#
