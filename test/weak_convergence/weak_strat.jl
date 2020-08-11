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
 Test Scalar SDEs (oop)
"""

@info "Scalar oop noise"

numtraj = Int(1e5) # in the paper they use 1e9
u₀ = 0.1
p = [1.5, 0.1]
f(u,p,t) = p[1]*u
g(u,p,t) = p[2]*u
dts = 1 .//2 .^(5:-1:0)
tspan = (0.0,1.0)

h1(z) = z

seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)



prob = SDEProblem(f,g,u₀,tspan, p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end]),false),
        prob_func = prob_func
        )
_solutions = @time generate_weak_solutions(ensemble_prob, RS1(), dts, numtraj, ensemblealg=EnsembleThreads())



errors = [LinearAlgebra.norm(Statistics.mean(sol.u) - u₀*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.3
println("RS1:", m)


seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, RS2(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u) - u₀*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.37

println("RS2:", m)

seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f,g,u₀,tspan, p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end]),false),
        prob_func = prob_func
        )


_solutions = @time generate_weak_solutions(ensemble_prob, NON(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.3

println("NON:", m)

numtraj = Int(1e5)
seed = 10
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, COM(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[2])/log(dts[end]/dts[2])
@test abs(m-2) < 0.3

println("COM:", m)


numtraj = Int(1e5)
seed = 10
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, NON2(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[2])/log(dts[end]/dts[2])
@test abs(m-2) < 0.3

println("NON2:", m)


"""
 Test Scalar SDEs (iip)
"""

@info "Scalar iip noise"

f!(du,u,p,t) = du[1] = p[1]*u[1]
g!(du,u,p,t) = du[1] = p[2]*u[1]

seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f!,g!,[u₀],tspan, p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end]),false),
        prob_func = prob_func
        )
_solutions = @time generate_weak_solutions(ensemble_prob, RS1(), dts,
  numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.3
println("RS1:", m)


seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f!,g!,[u₀],tspan, p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end]),false),
        prob_func = prob_func
        )
_solutions = @time generate_weak_solutions(ensemble_prob, RS2(), dts,
  numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.37

println("RS1:", m)

numtraj = Int(1e5)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f!,g!,[u₀],tspan, p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end]),false),
        prob_func = prob_func
        )
_solutions = @time generate_weak_solutions(ensemble_prob, NON(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.3

println("NON:", m)

numtraj = Int(1e5)
seed = 10
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, COM(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[2])/log(dts[end]/dts[2])
@test abs(m-2) < 0.3

println("COM:", m)

numtraj = Int(1e5)
seed = 10
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, NON2(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[2])/log(dts[end]/dts[2])
@test abs(m-2) < 0.3

println("NON2:", m)

"""
 Test Diagonal noise SDEs (iip)
"""

@info "Diagonal noise"


u₀ = [0.1,0.1]
function f3!(du,u,p,t)
  du[1] = 299//200*u[1]
  du[2] = 299//200*u[2]
end
function g3!(du,u,p,t)
  du[1] = 1//10*u[1]
  du[2] = 1//10*u[2]
end
dts = 1 .//2 .^(4:-1:0)
tspan = (0.0,1.0)

h3(z) = z^2 # == 1//10**exp(3//2*t) if h3(z) = z and  == 1//100**exp(301//100*t) if h3(z) = z^2 )

prob = SDEProblem(f3!,g3!,u₀,tspan)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h3(sol[end][1]),false),
        prob_func = prob_func
        )

numtraj = Int(4e4)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, RS1(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(301//100)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.37

println("RS1:", m)

numtraj = Int(4e4)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, RS2(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(301//100)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.3

println("RS2:", m)


numtraj = Int(1e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, NON(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(301//100)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.3

println("NON:", m)


numtraj = Int(6e4)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, COM(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(301//100)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.3

println("COM:", m)


numtraj = Int(6e4)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, NON2(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(301//100)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.3

println("NON2:", m)
