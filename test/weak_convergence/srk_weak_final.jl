"""
 Tests for https://arxiv.org/abs/1303.5103 with test problems as in the paper.
 DRI1 (and RI1)
"""


import Statistics # for mean values of trajectories
import LinearAlgebra # for the normn
using StochasticDiffEq
using Test
using Random
Random.seed!(100)
#using DiffEqGPU
#using Plots

"""
 Test Scalar SDEs (oop)
"""

numtraj = 5e6 # in the paper they use 1e9
u₀ = 0.0
f(u,p,t) = 1//2*u+sqrt(u^2+1)
g(u,p,t) = sqrt(u^2+1)
dts = 1 .//2 .^(4:-1:1)
tspan = (0.0,2.0)
prob = SDEProblem(f,g,u₀,tspan)

h1(z) = z^3-6*z^2+8*z

#analytical_sol(t) = E(f(X(t))) = E(h1(arsinh(X(t))) = t^3-3*t^2+2*t
#analytical_sol(2) = 0 and analytical_sol(1)=0
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(asinh(sol[end])),false)
        )
N = length(dts)
_solutions = @time [solve(ensemble_prob,
        DRI1();
        #ensemblealg = EnsembleGPUArray(); # EnsembleGPUArray or EnsembleDistributed()
        dt=dts[i],
        save_start=false,
        save_everystep=false,
        weak_timeseries_errors=false,
        weak_dense_errors=false,
        trajectories=numtraj) for i in 1:N]

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

#convergence_plot = plot(dts, errors, xaxis=:log, yaxis=:log)
#savefig(convergence_plot, "DRI1-CPU-final-"*string(numtraj)*".pdf")
println("DRI1:", m)



_solutions = @time [solve(ensemble_prob,
        RI1();
        #ensemblealg = EnsembleGPUArray(); # EnsembleGPUArray or EnsembleDistributed()
        dt=dts[i],
        save_start=false,
        save_everystep=false,
        weak_timeseries_errors=false,
        weak_dense_errors=false,
        trajectories=numtraj) for i in 1:N]

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

# convergence_plot = plot(dts, errors, xaxis=:log, yaxis=:log)
# savefig(convergence_plot, "RI1-CPU-final-"*string(numtraj)*".pdf")
println("RI1:", m)


"""
 Test Scalar SDEs (iip)
"""

numtraj = 5e6
u₀ = [0.0]
f1!(du,u,p,t) = @.(du = 1//2*u+sqrt(u^2 +1))
g1!(du,u,p,t) = @.(du = sqrt(u^2 +1))
dts = 1 .//2 .^(4:-1:1)
tspan = (0.0,2.0)
prob = SDEProblem(f1!,g1!,u₀,tspan)

h1(z) = z^3-6*z^2+8*z

ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(asinh(sol[end][1])),false)
        )
N = length(dts)
_solutions = @time [solve(ensemble_prob,
        DRI1();
        dt=dts[i],
        save_start=false,
        save_everystep=false,
        weak_timeseries_errors=false,
        weak_dense_errors=false,
        trajectories=numtraj) for i in 1:N]

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

# convergence_plot = plot(dts, errors, xaxis=:log, yaxis=:log)
# savefig(convergence_plot, "Scalar-DRI1-CPU-final-"*string(numtraj)*".pdf")
println("DRI1:", m)

_solutions = @time [solve(ensemble_prob,
        RI1();
        dt=dts[i],
        save_start=false,
        save_everystep=false,
        weak_timeseries_errors=false,
        weak_dense_errors=false,
        trajectories=numtraj) for i in 1:N]

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("RI1:", m)

"""
 Test non-commutative noise SDEs (iip)
"""

numtraj = 1e6
u₀ = [1.0,1.0]
function f2!(du,u,p,t)
  du[1] = -273//512*u[1]
  du[2] = -1//160*u[1]-(-785//512+sqrt(2)/8)*u[2]
end
function g2!(du,u,p,t)
  du[1,1] = 1//4*u[1]
  du[1,2] = 1//16*u[1]
  du[2,1] = (1-2*sqrt(2))/4*u[1]
  du[2,2] = 1//10*u[1]+1//16*u[2]
end
dts = 1 .//2 .^(3:-1:0)
tspan = (0.0,10.0)
prob = SDEProblem(f2!,g2!,u₀,tspan,noise_rate_prototype=zeros(2,2))

h2(z) = z^2 # but apply it only to u[1]

ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h2(sol[end][1]),false)
        )
N = length(dts)
_solutions = @time [solve(ensemble_prob,
        DRI1();
        dt=dts[i],
        save_start=false,
        save_everystep=false,
        weak_timeseries_errors=false,
        weak_dense_errors=false,
        trajectories=numtraj) for i in 1:N]

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-exp(-10.0)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

# convergence_plot = plot(dts, errors, xaxis=:log, yaxis=:log)
# savefig(convergence_plot, "ND-DRI1-CPU-final-"*string(numtraj)*".pdf")
println("DRI1:", m)

_solutions = @time [solve(ensemble_prob,
        RI1();
        dt=dts[i],
        save_start=false,
        save_everystep=false,
        weak_timeseries_errors=false,
        weak_dense_errors=false,
        trajectories=numtraj) for i in 1:N]

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-exp(-10.0)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("RI1:", m)

"""
 Test Diagonal noise SDEs (iip), SIAM Journal on Numerical Analysis, 47 (2009), pp. 1713–1738
"""

numtraj = 1e6
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
prob = SDEProblem(f3!,g3!,u₀,tspan)

h3(z) = z^2 # == 1//10**exp(3//2*t) if h3(z) = z and  == 1//100**exp(301//100*t) if h3(z) = z^2 )

ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h3(sol[end][1]),false)
        )
N = length(dts)
_solutions = @time [solve(ensemble_prob,
        DRI1();
        dt=dts[i],
        save_start=false,
        save_everystep=false,
        weak_timeseries_errors=false,
        weak_dense_errors=false,
        trajectories=numtraj) for i in 1:N]

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(301//100)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.5

println("DRI1:", m)

#convergence_plot = plot(dts, errors, xaxis=:log, yaxis=:log)
#savefig(convergence_plot, "RI1-CPU-final-"*string(numtraj)*".pdf")
#println(m)

_solutions = @time [solve(ensemble_prob,
        RI1();
        dt=dts[i],
        save_start=false,
        save_everystep=false,
        weak_timeseries_errors=false,
        weak_dense_errors=false,
        trajectories=numtraj) for i in 1:N]

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(301//100)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.5
println("RI1:", m)
