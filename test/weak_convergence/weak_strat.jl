"""
 Tests for  https://link.springer.com/article/10.1007/s10543-007-0130-3 with test problems as in the paper.
 RS1, RS2
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
u₀ = 0.5
p = [1, 1//2]
f(u,p,t) = p[1]*u
g(u,p,t) = p[2]*u
dts = 1 .//2 .^(6:-1:1)
tspan = (0.0,1.0) # 2.0 in paper


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
@test -(m-2) < 0.3

using Plots
convergence_plot = plot(dts, errors, xaxis=:log, yaxis=:log)
savefig(convergence_plot, "RS1-final-"*string(numtraj)*".pdf")

f_strat(u0,p,t,W) = @.(u0*exp(1.5t+0.1W))
prob_test = SDEProblem(SDEFunction(f,g,analytic=f_strat),g,0.1,(0.0,1.0))

ensemble_test = EnsembleProblem(prob_test;
        #output_func = (sol,i) -> (h1(sol[end]),false),
        prob_func = prob_func
        )

sol_test = solve(ensemble_test, EulerHeun(), dt=dts[1], trajectories=Int(1e5))

summ = EnsembleSummary(sol_test,0:0.01:1)

plot(summ,fillalpha=0.5)


analy = Statistics.mean([ sol.u_analytic for sol in sol_test])
ts =  Statistics.mean([ sol.t for sol in sol_test])

exp_analy = @. u₀*exp(1.5*(ts))

plot(ts, analy)
plot!(ts, exp_analy)

plot(ts, analy-exp_analy)

exp2 = @. u₀*exp((1.5+0.5*0.1^2)*(ts))
plot(ts, analy-exp2)


solve(prob, RS2(), dt=1//4)

f(0.5,1.0,0.25)

println("RS1:", m)



numtraj = Int(1e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, RS2(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [abs(LinearAlgebra.norm(Statistics.mean(sol.u)) - u₀*exp(1.0*1.5)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("RS2:", m)

"""
 Test Scalar SDEs (iip)
"""

@info "Scalar iip noise"


"""
 Test non-commutative noise SDEs (iip)
"""

@info "Non-commutative noise"




"""
 Test Diagonal noise SDEs (iip), SIAM Journal on Numerical Analysis, 47 (2009), pp. 1713–1738
"""

@info "Diagonal noise"
