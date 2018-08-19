using StochasticDiffEq, DiffEqProblemLibrary, DiffEqDevTools, Test, LinearAlgebra, Random

Random.seed!(100)
const σ_const = 0.87
const μ = 1.01

u0 = rand(2)
A = [-2.0 1.0;1.0 -2.0]
B = Diagonal([σ_const for i in 1:2])

function f_nondiag(u,p,t)
  return A*u + 1.01*u
end

function f_nondiag(du,u,p,t)
  mul!(du,A,u)
  du .+= 1.01u
end

function f_nondiag(::Type{Val{:analytic}},u0,p,t,W)
 tmp = (A+1.01I-(B^2))*t + B*sum(W)
 exp(tmp)*u0
end

function g_nondiag(u,p,t)
  du = zeros(2,2)
  du[1,1] = σ_const*u[1]
  du[1,2] = σ_const*u[1]
  du[2,1] = σ_const*u[2]
  du[2,2] = σ_const*u[2]
  return du
end

function g_nondiag(du,u,p,t)
  du[1,1] = σ_const*u[1]
  du[1,2] = σ_const*u[1]
  du[2,1] = σ_const*u[2]
  du[2,2] = σ_const*u[2]
end

coeff = 2*σ_const^2 #To not compute the same coefficient in the sde.
function ggprime(u, p, t)
  return coeff*u
end

function ggprime(du, u, p, t)
  du .= coeff*u
end

prob = SDEProblem{false}(f_nondiag,g_nondiag,u0,(0.0,1.0),noise_rate_prototype=zeros(2,2))
probiip = SDEProblem{true}(f_nondiag,g_nondiag,u0,(0.0,1.0),noise_rate_prototype=zeros(2,2))

## Just solve to test compatibility
IEM = solve(probiip,ImplicitEM())
IEM = solve(probiip,ImplicitEulerHeun())
IEM = solve(probiip,ISSEM())
IEM = solve(probiip,ISSEulerHeun())

## use BenchmarkTools to benchmark speed here if desired.
dt = 0.002
seed = 1
# solEMiip = solve(probiip, EM(), dt=dt, seed=seed)
# solRKMil = solve(probiip, RKMilCommute(), dt=dt, seed=seed)
solPCEuler = solve(prob, PCEuler(ggprime), dt=dt, seed=seed)
solPCEuleriip = solve(probiip, PCEuler(ggprime), dt=dt, seed=seed)

@test solPCEuler.u ≈ solPCEuleriip.u atol=1e-10

## plot to see performance of PCEuler
# first_elem(x) = x[1]
# using PyPlot
# figure()
# plot(solEMiip.t, first_elem.(solEMiip.u_analytic), label="True", "--")
# plot(solEMiip.t, first_elem.(solEMiip.u), label="EM")
# plot(solRKMil.t, first_elem.(solRKMil.u), label="RKMil")
# plot(solPCEuleriip.t, first_elem.(solPCEuleriip.u), label="PC")
# legend()

##
dts = (1/2) .^ (10:-1:5) #14->7 good plot
numMonte = 50
simEM = test_convergence(dts,probiip,EM(),numMonte=numMonte)
simPCEuler = test_convergence(dts,probiip,PCEuler(ggprime),numMonte=numMonte)
#simRKMil = test_convergence(dts,probiip,RKMilCommute(),numMonte=numMonte)
@test all(simPCEuler.errors[:l2] .< simEM.errors[:l2])

## Plotting script to see the order 1 scaling of PCEuler
# figure()
# loglog(dts, simEM.errors[:l2], label="EM")
# loglog(dts, simRKMil.errors[:l2], label="RKMil")
# loglog(dts, simPCEuler.errors[:l2], label="PC")
# ylabel("l2 error")
# xlabel("dt")
# legend()

f_morenoise = (du,u,p,t) -> du.=1.01u
g_morenoise = function (du,u,p,t)
  du[1,1] = 0.3u[1]
  du[1,2] = 0.6u[1]
  du[1,3] = 0.9u[1]
  du[1,4] = 0.12u[2]
  du[2,1] = 1.2u[1]
  du[2,2] = 0.2u[2]
  du[2,3] = 0.3u[2]
  du[2,4] = 1.8u[2]
end
prob = SDEProblem(f_morenoise,g_morenoise,ones(2),(0.0,1.0),
                  noise_rate_prototype=zeros(2,4))

sol =solve(prob,dt=1/2^(3),EM())
sol =solve(prob,dt=1/2^(3),ISSEM())
sol =solve(prob,dt=1/2^(3),ImplicitEM())
sol =solve(prob,dt=1/2^(3),EulerHeun())
sol =solve(prob,dt=1/2^(3),ImplicitEulerHeun())
sol =solve(prob,dt=1/2^(3),ISSEulerHeun())
