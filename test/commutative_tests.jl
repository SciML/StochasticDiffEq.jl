using StochasticDiffEq, DiffEqProblemLibrary, DiffEqDevTools, Base.Test
srand(100)
dts = 1./2.^(10:-1:2) #14->7 good plot

using SpecialMatrices
const σ_const = 0.87
const μ = 1.01

u0 = rand(2)
A = full(Strang(2))
B = Diagonal([σ_const for i in 1:2])

function f_commute(t,u,du)
  A_mul_B!(du,A,u)
  du .+= 1.01u
end
function (p::typeof(f_commute))(::Type{Val{:analytic}},t,u0,W)
 tmp = (A+1.01I-(B^2))*t + B*sum(W)
 expm(tmp)*u0
end
function σ(t,u,du)
  du[1,1] = σ_const*u[1]
  du[1,2] = σ_const*u[1]
  du[2,1] = σ_const*u[2]
  du[2,2] = σ_const*u[2]
end

prob = SDEProblem(f_commute,σ,u0,(0.0,1.0),noise_rate_prototype=rand(2,2))

sol = solve(prob,RKMilCommute(),dt=1/2^(8))
sol = solve(prob,EM(),dt=1/2^(10))

dts = 1./2.^(10:-1:3) #14->7 good plot
sim2 = test_convergence(dts,prob,EM(),numMonte=Int(1e2))
sim2 = test_convergence(dts,prob,RKMilCommute(),numMonte=Int(1e2))
