"""
 Tests for  https://link.springer.com/article/10.1007/s10543-007-0130-3 with test problems as in the paper.
 RS1, RS2
 and for https://www.sciencedirect.com/science/article/pii/S0377042706003906
 NON, COM, NON2
"""


import Statistics # for mean values of trajectories
import LinearAlgebra # for the normn
using StochasticDiffEq
using Test
using Random
#using DiffEqGPU

function generate_weak_solutions(prob, alg, dts, numtraj; ensemblealg=EnsembleThreads())
  sols = []
  for i in 1:length(dts)
    sol = solve(prob,alg;ensemblealg=ensemblealg,dt=dts[i],adaptive=false,save_start=false,save_everystep=false,weak_timeseries_errors=false,weak_dense_errors=false,trajectories=Int(numtraj))
    println(i)
    push!(sols,sol)
  end
  return sols
end

function prob_func(prob, i, repeat)
    remake(prob,seed=seeds[i])
end

"""
 Test non-commutative noise SDEs (iip)
"""

@info "Non-commutative noise"

u₀ = [0.1,0.1]
function f2!(du,u,p,t)
  du[1] = 5//4*u[2]-5//4*u[1]
  du[2] = 1//4*u[1]-1//4*u[2]
end
function g2!(du,u,p,t)
  du[1,1] = sqrt(3)/2*(u[1]-u[2])
  du[1,2] = 1//2*(u[1]+u[2])
  #du[2,1] = 0
  du[2,2] = u[1]
end
dts = 1 .//2 .^(4:-1:1)
tspan = (0.0,1.0)

h2(z) = z^2 # E(x_i) = 1/10 exp(1/2t) or E(x_1* x_2) = 1/100 exp(2t)

prob = SDEProblem(f2!,g2!,u₀,tspan,noise_rate_prototype=zeros(2,2))
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h2(sol[end][1]),false),
        prob_func = prob_func
        )

numtraj = Int(1e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, RS1(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(2)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.3

println("RS1:", m)


numtraj = Int(1e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, RS2(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(2)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.3

println("RS2:", m)

numtraj = Int(1e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, NON(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(2)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.3

println("NON:", m)

numtraj = Int(5e5)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, COM(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(2)) for sol in _solutions]
m = log(errors[end]/errors[2])/log(dts[end]/dts[2])
@test abs(m-2) < 0.3 # tests are passing; problem might be not hard enough..

println("COM:", m)


numtraj = Int(5e5)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, NON2(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(2)) for sol in _solutions]
m = log(errors[end]/errors[2])/log(dts[end]/dts[2])
@test abs(m-2) < 0.3 # tests are passing; problem might be not hard enough..

println("NON2:", m)
