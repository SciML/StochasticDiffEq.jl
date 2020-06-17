using StochasticDiffEq, Test, Random
using DiffEqGPU
using CuArrays
import Statistics # for mean values of trajectories
import LinearAlgebra # for the normn

function weak_error(prob, alg, numtraj, expected_value;abstol=1,reltol=0,ensemblealg=EnsembleCPUArray())
  sol = solve(prob,alg;ensemblealg=ensemblealg,dt=0.0625f0, abstol=abstol,reltol=reltol,
      save_start=false,save_everystep=false,weak_timeseries_errors=false,weak_dense_errors=false,
      trajectories=Int(numtraj))
  return LinearAlgebra.norm(Statistics.mean(sol.u)-expected_value)
end

function prob_func(prob, i, repeat)
    remake(prob,seed=seeds[i])
end

# prob 1
u₀ = [0.0f0]
tspan = (0.0f0,2.0f0)
h1(z) = z^3-6*z^2+8*z

function f1!(du,u,p,t)
 @inbounds begin
     du[1] = 1//2*u[1]+sqrt(u[1]^2 +1)
 end
 nothing
end

function g1!(du,u,p,t)
 @inbounds begin
     du[1] = sqrt(u[1]^2 +1)
 end
 nothing
end

prob1 = SDEProblem(f1!,g1!,u₀,tspan)
ensemble_prob1 = EnsembleProblem(prob1;
        output_func = (sol,i) -> (h1(asinh(sol[end][1])),false),
        prob_func = prob_func
        )

probs = Vector{EnsembleProblem}(undef, 1)
probs[1] = ensemble_prob1


numtraj = Int(1e5)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)



for i in 1:1
  err1 = weak_error(probs[i],DRI1(),numtraj,0.0f0,abstol=1,reltol=0)
  @show err1
  err2 = weak_error(probs[i],DRI1(),numtraj,0.0f0,abstol=1e-1,reltol=0)
  @show err2
  err3 = weak_error(probs[i],DRI1(),numtraj,0.0f0,abstol=1e-2,reltol=0)
  @show err3
  @test err1 > err2
  @test err2 > err3
end




tsave = 0.0f0:0.05f0:2f0
ensemble_test = EnsembleProblem(prob1;
        prob_func = prob_func
        )
sol = solve(ensemble_test,DRI1();ensemblealg=EnsembleCPUArray(),dt=0.001f0,#adaptive=false,
    abstol=1e-2,reltol=0,trajectories=Int(numtraj),saveat=tsave)



fsol = h1.(asinh.(sol))
computed_exp = Statistics.mean(fsol, dims=3)[1,:,1]
var = Statistics.std(fsol, dims=3)[1,:,1]/sqrt(numtraj)
true_exp = tsave.^3-3*tsave.^2+2*tsave
#
# using Plots; plt=plot(Array(tsave),true_exp, label="true", xaxis="time t",yaxis="E(X(t))", title="dt=0.001 adaptive", lw=3.0)
# plot!(Array(tsave), computed_exp, ribbon = var, label="DRI1")
# #
# savefig(plt, "adaptive-1e-2.png")
