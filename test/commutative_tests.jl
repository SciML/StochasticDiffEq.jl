@everywhere using StochasticDiffEq, DiffEqProblemLibrary, DiffEqDevTools, Base.Test
srand(100)
dts = 1./2.^(10:-1:2) #14->7 good plot

using SpecialMatrices
const σ_const = 0.87
const μ = 1.01

u0 = rand(2)
A = Strang(2)
B = Diagonal([σ_const for i in 1:2])

function f(t,u,du)
  A_mul_B!(du,A,u)
  du .+= 1.01u
end
function (p::typeof(f))(::Type{Val{:analytic}},t,u0,W)
 tmp = (A+1.01I-(B^2))*t + B*sum(W)
 expm(tmp)*u0
end
function σ(t,u,du)
  du[1,1] = σ_const*u[1]
  du[1,2] = σ_const*u[1]
  du[2,1] = σ_const*u[2]
  du[2,2] = σ_const*u[2]
end

prob = SDEProblem(f,σ,u0,(0.0,1.0),noise_rate_prototype=rand(2,2))

sol = solve(prob,RKMilCommute(),dt=1/2^(8))
sol = solve(prob,EM(),dt=1/2^(10))

dts = 1./2.^(10:-1:3) #14->7 good plot
sim2 = test_convergence(dts,prob,EM(),numMonte=Int(1e2))
sim2 = test_convergence(dts,prob,RKMilCommute(),numMonte=Int(4e2))

f1 = (t,u,du) -> A
function (p::typeof(f1))(::Type{Val{:analytic}},t,u0,W)
 tmp = (A+1.01I-(B^2))*t + B*sum(W)
 expm(tmp)*u0
end
f2 = (t,u,du) -> du .= μ .* u
function σ(t,u,du)
  du[1,1] = σ_const*u[1]
  du[1,2] = σ_const*u[1]
  du[2,1] = σ_const*u[2]
  du[2,2] = σ_const*u[2]
end
prob = SDEProblem((f1,f2),σ,u0,(0.0,1.0),noise_rate_prototype=rand(2,2))
sol = solve(prob,IIF1M(),dt=1/10)
sol = solve(prob,IIF1Mil(),dt=1/10)

sim  = test_convergence(dts,prob,IIF1M(),numMonte=Int(5e1))
@test abs(sim.𝒪est[:l2]-0.5) < 0.2

sim  = test_convergence(dts,prob,IIF1Mil(),numMonte=Int(4e2))
@test abs(sim.𝒪est[:l2]-1.0) < 0.2

u0 = rand(5)
A = Strang(5)
B = Diagonal([σ_const for i in 1:5])
function f(t,u,du)
  A_mul_B!(du,A,u)
end
function σ(t,u,du)
  du[1,1] = σ_const*u[1]
  du[1,2] = σ_const*u[1]
  du[1,3] = σ_const*u[1]
  du[1,4] = σ_const*u[1]
  du[1,5] = σ_const*u[1]
  du[2,1] = σ_const*u[2]
  du[2,2] = σ_const*u[2]
  du[2,3] = σ_const*u[2]
  du[2,4] = σ_const*u[2]
  du[2,5] = σ_const*u[2]
  du[3,1] = σ_const*u[3]
  du[3,2] = σ_const*u[3]
  du[3,3] = σ_const*u[3]
  du[3,4] = σ_const*u[3]
  du[3,5] = σ_const*u[3]
  du[4,1] = σ_const*u[4]
  du[4,2] = σ_const*u[4]
  du[4,3] = σ_const*u[4]
  du[4,4] = σ_const*u[4]
  du[4,5] = σ_const*u[4]
  du[5,1] = σ_const*u[5]
  du[5,2] = σ_const*u[5]
  du[5,3] = σ_const*u[5]
  du[5,4] = σ_const*u[5]
  du[5,5] = σ_const*u[5]
end
function (p::typeof(f))(::Type{Val{:analytic}},t,u0,W)
 tmp = (A+1.01I-(5/2)*(B^2))*t + B*sum(W)
 expm(tmp)*u0
end
A*B == B*A

prob = SDEProblem(f,σ,u0,(0.0,1.0),noise_rate_prototype=rand(5,5))
sol = solve(prob,RKMilCommute(),dt=1/2^(8))
sol = solve(prob,EM(),dt=1/2^(10))

dts = 1./2.^(10:-1:3) #14->7 good plot
sim2 = test_convergence(dts,prob,EM(),numMonte=Int(1e2))
sim2 = test_convergence(dts,prob,RKMilCommute(),numMonte=Int(1e2))


f1 = (t,u,du) -> A
function (p::typeof(f1))(::Type{Val{:analytic}},t,u0,W)
 tmp = (A+1.01I-(5/2)*(B^2))*t + B*sum(W)
 expm(tmp)*u0
end
f2 = (t,u,du) -> du .= μ .* u
prob = SDEProblem((f1,f2),σ,u0,(0.0,1.0),noise_rate_prototype=rand(5,5))
sol = solve(prob,IIF1M(),dt=1/10)
sol = solve(prob,IIF1Mil(),dt=1/10)

dts = 1./2.^(10:-1:3) #14->7 good plot
sim2 = test_convergence(dts,prob,IIF1M(),numMonte=Int(1e2))
sim2 = test_convergence(dts,prob,IIF1Mil(),numMonte=Int(1e2))

################################################################################

u0 = rand(5)
A = Strang(5)
B = Diagonal([σ_const for i in 1:5])
function f(t,u,du)
  A_mul_B!(du,A,u)
end
function σ(t,u,du)
  du[1,1] = σ_const*u[1]
  du[1,2] = σ_const*u[1]
  du[1,3] = σ_const*u[1]
  du[1,4] = σ_const*u[1]
  du[1,5] = σ_const*u[1]
  du[2,1] = σ_const*u[2]
  du[2,2] = σ_const*u[2]
  du[2,3] = σ_const*u[2]
  du[2,4] = σ_const*u[2]
  du[2,5] = σ_const*u[2]
  du[3,1] = σ_const*u[3]
  du[3,2] = σ_const*u[3]
  du[3,3] = σ_const*u[3]
  du[3,4] = σ_const*u[3]
  du[3,5] = σ_const*u[3]
  du[4,1] = σ_const*u[4]
  du[4,2] = σ_const*u[4]
  du[4,3] = σ_const*u[4]
  du[4,4] = σ_const*u[4]
  du[4,5] = σ_const*u[4]
  du[5,1] = σ_const*u[5]
  du[5,2] = σ_const*u[5]
  du[5,3] = σ_const*u[5]
  du[5,4] = σ_const*u[5]
  du[5,5] = σ_const*u[5]
end
prob = SDEProblem(f,σ,u0,(0.0,1.0),noise_rate_prototype=rand(5,5))

using DiffEqNoiseProcess
dts = 1./2.^(6:-1:3) #14->7 good plot
N = length(dts)
_solutions = []
numMonte = 100
adaptive = false
alg = EM()
test_dt = 1./2.^12
for i in 1:N
  tmp_solutions = []
  for j in 1:numMonte
    sol = solve(prob,alg;dt=dts[i],adaptive=adaptive);
    prob2 = SDEProblem(prob.f,prob.g,prob.u0,prob.tspan,
                       noise=NoiseWrapper(sol.W),
                       noise_rate_prototype=prob.noise_rate_prototype);
    true_sol = solve(prob2,alg;adaptive=adaptive,dt=test_dt);
    err_sol = appxtrue(sol,true_sol)
    push!(tmp_solutions,err_sol)
  end
  push!(_solutions,MonteCarloSolution(tmp_solutions,0.0,true))
end



test_dt = 1./2.^(12)
dts = 1./2.^(10:-1:3) #14->7 good plot
dt = 0.001
t = 0:test_dt:1
brownian_values = cumsum([[zeros(5)];[sqrt(test_dt)*randn(5) for i in 1:length(t)-1]])
W = NoiseGrid(t,brownian_values)
W(0.9)

u0 = rand(5)
function f(t,u,du)
  du.=u
end
function σ(t,u,du)
  du.=1.0u
end

W = NoiseGrid(t,brownian_values)
prob = SDEProblem(f,σ,u0,(0.0,1.0),noise=W)
sol1 = solve(prob,EM(),dt=1/10)
W = NoiseGrid(t,brownian_values)
prob = SDEProblem(f,σ,u0,(0.0,1.0),noise=W)
sol2 = solve(prob,EM(),dt=1/100)
using Plots; plot(sol1); plot!(sol2)

sim1 = DiffEqDevTools.analyticless_test_convergence(dts,prob,EM(),test_dt,numMonte=Int(1e2))
sim2 = DiffEqDevTools.analyticless_test_convergence(dts,prob,RKMil(),test_dt,numMonte=Int(1e2))













u0 = rand(5)
A = Strang(5)
B = Diagonal([σ_const for i in 1:5])
function f(t,u,du)
  A_mul_B!(du,A,u)
end
function σ(t,u,du)
  du[1,1] = σ_const*u[1]
  du[1,2] = σ_const*u[1]
  du[1,3] = σ_const*u[1]
  du[1,4] = σ_const*u[1]
  du[1,5] = σ_const*u[1]
  du[2,1] = σ_const*u[2]
  du[2,2] = σ_const*u[2]
  du[2,3] = σ_const*u[2]
  du[2,4] = σ_const*u[2]
  du[2,5] = σ_const*u[2]
  du[3,1] = σ_const*u[3]
  du[3,2] = σ_const*u[3]
  du[3,3] = σ_const*u[3]
  du[3,4] = σ_const*u[3]
  du[3,5] = σ_const*u[3]
  du[4,1] = σ_const*u[4]
  du[4,2] = σ_const*u[4]
  du[4,3] = σ_const*u[4]
  du[4,4] = σ_const*u[4]
  du[4,5] = σ_const*u[4]
  du[5,1] = σ_const*u[5]
  du[5,2] = σ_const*u[5]
  du[5,3] = σ_const*u[5]
  du[5,4] = σ_const*u[5]
  du[5,5] = σ_const*u[5]
end
prob = SDEProblem(f,σ,u0,(0.0,1.0),noise_rate_prototype=rand(5,5))

sim3 = DiffEqDevTools.analyticless_test_convergence(dts,prob,EM(),test_dt,numMonte=Int(1e2))
sim4 = DiffEqDevTools.analyticless_test_convergence(dts,prob,RKMilCommute(),test_dt,numMonte=Int(1e2))

f1 = (t,u,du) -> A
f2 = (t,u,du) -> du .= μ .* u
prob = SDEProblem((f1,f2),σ,u0,(0.0,1.0),noise_rate_prototype=rand(5,5))

sol = solve(prob,IIF1M(),dt=1/10)
sol = solve(prob,IIF1Mil(),dt=1/10)

sim5 = DiffEqDevTools.analyticless_test_convergence(dts,prob,IIF1M(),test_dt,numMonte=Int(1e2))



















sim6 = DiffEqDevTools.analyticless_test_convergence(dts,prob,IIF1Mil(),test_dt,numMonte=Int(1e2))
