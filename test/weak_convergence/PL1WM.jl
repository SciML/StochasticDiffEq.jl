"""
 Tests for PL1WM
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
    sol = solve(prob,alg;ensemblealg=ensemblealg,dt=dts[i],save_start=false,save_everystep=false,weak_timeseries_errors=false,weak_dense_errors=false,trajectories=Int(numtraj))
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

# PC exercise 14.2.2
numtraj = Int(1e4)
u₀ = 0.1
f(u,p,t) = p[1]*u
g(u,p,t) = p[2]*u
dts = 1 .//2 .^(6:-1:3)
tspan = (0.0,1.0)
p = [3//2,1//100]

h1(z) = z
#analytical_sol(t) = E(f(X(t))) =

seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f,g,u₀,tspan,p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end]),false),
        prob_func = prob_func
        )
_solutions = @time generate_weak_solutions(ensemble_prob, PL1WM(), dts, numtraj, ensemblealg=EnsembleThreads())

errors =  [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

#using Plots; convergence_plot = plot(dts, errors, xaxis=:log, yaxis=:log)
#savefig(convergence_plot, "PL1WM-scalar.png")
println("PL1WM:", m)


"""
 Test Scalar SDEs (iip)
"""

@info "Scalar iip noise"

u₀ = [0.1]
f1!(du,u,p,t) =  (du[1] = p[1]*u[1])
g1!(du,u,p,t) =  (du[1] = p[2]*u[1])
dts = 1 .//2 .^(6:-1:3)
tspan = (0.0,1.0)
p = [3//2,1//100]


prob = SDEProblem(f1!,g1!,u₀,tspan,p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end][1]),false),
        prob_func = prob_func
        )


numtraj = Int(1e4)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)


_solutions = @time generate_weak_solutions(ensemble_prob, PL1WM(), dts, numtraj, ensemblealg=EnsembleThreads())

errors =  [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("PL1WM:", m)

"""
 Test non-commutative noise SDEs (iip)
"""

@info "Non-commutative noise"

u₀ = [1/8,1/8,1/1,1/8]
function f2!(du,u,p,t)
  du[1] = 243//154*u[1] - 27//77*u[2] + 23//154*u[3] - 65//154*u[4]
  du[2] = 27//77*u[1] - 243//154*u[2] + 65//154*u[3] - 23//154*u[4]
  du[3] = 5//154*u[1] - 61//154*u[2] + 162//77*u[3] - 36//77*u[4]
  du[4] = 61//154*u[1] - 5//154*u[2] + 36//77*u[3] - 162//77*u[4]
end
function g2!(du,u,p,t)
  du[1,1] = 1//9*sqrt(u[2]^2+u[3]^2+2//23)*1//13
  du[1,2] = 1//8*sqrt(u[4]^2+u[1]^2+1//11)*1//14
  du[1,3] = p[1]*1//12*sqrt(u[1]^2+u[2]^2+1//9)*1//6
  du[1,4] = p[1]*1//14*sqrt(u[3]^2+u[4]^2+3//29)*1//8
  du[2,1] = 1//9*sqrt(u[2]^2+u[3]^2+2//23)*1//14
  du[2,2] = 1//8*sqrt(u[4]^2+u[1]^2+1//11)*1//16
  du[2,3] = p[1]*1//12*sqrt(u[1]^2+u[2]^2+1//9)*1//5
  du[2,4] = p[1]*1//14*sqrt(u[3]^2+u[4]^2+3//29)*1//9
  du[3,1] = 1//9*sqrt(u[2]^2+u[3]^2+2//23)*1//13
  du[3,2] = 1//8*sqrt(u[4]^2+u[1]^2+1//11)*1//16
  du[3,3] = p[1]*1//12*sqrt(u[1]^2+u[2]^2+1//9)*1//5
  du[3,4] = p[1]*1//14*sqrt(u[3]^2+u[4]^2+3//29)*1//8
  du[4,1] = 1//9*sqrt(u[2]^2+u[3]^2+2//23)*1//15
  du[4,2] = 1//8*sqrt(u[4]^2+u[1]^2+1//11)*1//12
  du[4,3] = p[1]*1//12*sqrt(u[1]^2+u[2]^2+1//9)*1//6
  du[4,4] = p[1]*1//14*sqrt(u[3]^2+u[4]^2+3//29)*1//9
end
dts = 1 .//2 .^(6:-1:1)
tspan = (0.0,1.0)
p = [1]
h2(z) = z
# solution: E(X^i) = 1/8 exp(2*T), for i=1,2,4; E(X^3) = exp(2*T)

prob = SDEProblem(f2!,g2!,u₀,tspan,p,noise_rate_prototype=zeros(4,4))
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h2(sol[end][1]),false),
        prob_func = prob_func
        )

numtraj = Int(1e4)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)


_solutions = @time generate_weak_solutions(ensemble_prob, PL1WM(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-u₀[1]*exp(2*1.0)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.33

println("PL1WM:", m)



"""
 Test Diagonal noise SDEs (iip), SIAM Journal on Numerical Analysis, 47 (2009), pp. 1713–1738
"""

@info "Diagonal noise"

u₀ = [0.1,0.1]
function f3!(du,u,p,t)
  du[1] = 3//2*u[1]
  du[2] = 3//2*u[2]
end
function g3!(du,u,p,t)
  du[1] = 1//10*u[1]
  du[2] = 1//10*u[1]
end
dts = 1 .//2 .^(3:-1:0)
tspan = (0.0,1.0)

h3(z) = z^2 # == 1//10**exp(3//2*t) if h3(z) = z and  == 1//100**exp(301//100*t) if h3(z) = z^2 )

prob = SDEProblem(f3!,g3!,u₀,tspan)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h3(sol[end][1]),false),
        prob_func = prob_func
        )

numtraj = Int(5e4)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, PL1WM(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(301//100)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("PL1WM:", m)



"""
 Test Additive noise SDEs Kloeden & Platen Exercise 14.4.1
"""

numtraj = Int(1e3)
u₀ = 0.1
fadd(u,p,t) = p[1]*u
gadd(u,p,t) = p[2]
dts = 1 .//2 .^(4:-1:0)
tspan = (0.0,1.0)
p = [2//2,1//100]

seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(fadd,gadd,u₀,tspan,p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (sol[end],false),
        prob_func = prob_func
        )
_solutions = @time generate_weak_solutions(ensemble_prob, PL1WM(), dts, numtraj, ensemblealg=EnsembleThreads())
_solutions1 = @time generate_weak_solutions(ensemble_prob, PL1WMA(), dts, numtraj, ensemblealg=EnsembleThreads())

errors =  [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3
@test minimum(_solutions .≈ _solutions1)

println("PL1WM:", m)

#inplace

u₀ = [0.1]
fadd!(du,u,p,t) = @.(du = p[1]*u)
gadd!(du,u,p,t) = @.(du = p[2])


prob = SDEProblem(fadd!,gadd!,u₀,tspan,p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (sol[end],false),
        prob_func = prob_func
        )
_solutions = @time generate_weak_solutions(ensemble_prob, PL1WM(), dts, numtraj, ensemblealg=EnsembleThreads())
_solutions1 = @time generate_weak_solutions(ensemble_prob, PL1WMA(), dts, numtraj, ensemblealg=EnsembleThreads())

errors =  [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3
@test minimum(_solutions .≈ _solutions1)
