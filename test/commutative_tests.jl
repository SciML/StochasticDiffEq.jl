@everywhere using StochasticDiffEq, DiffEqProblemLibrary, DiffEqDevTools, Base.Test
srand(100)
dts = 1./2.^(10:-1:2) #14->7 good plot

using SpecialMatrices
const σ_const = 0.87
const μ = 1.01

u0 = rand(2)
A = Strang(2)
B = [σ_const 0
    0 σ_const]

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

A = Strang(5)
f1 = (t,u,du) -> A
function (p::typeof(f1))(::Type{Val{:analytic}},t,u0,W)
 tmp = (A+1.01I-(B^2))*t + B*sum(W)
 expm(tmp)*u0
end
f2 = (t,u,du) -> du .= μ .* u
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
B = Diagonal([σ_const for i in 1:5])
u0 = rand(5)
prob = SDEProblem((f1,f2),σ,u0,(0.0,1.0),noise_rate_prototype=rand(5,5))
sol = solve(prob,IIF1M(),dt=1/10)
sol = solve(prob,IIF1Mil(),dt=1/10)

sim  = test_convergence(dts,prob,IIF1M(),numMonte=Int(5e1))
@test abs(sim.𝒪est[:l2]-0.5) < 0.2

sim  = test_convergence(dts,prob,IIF1Mil(),numMonte=Int(4e2))
@test abs(sim.𝒪est[:l2]-1.0) < 0.2
