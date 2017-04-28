using DiffEqBase, StochasticDiffEq, DiffEqNoiseProcess, Base.Test, DiffEqDevTools, SpecialMatrices
const μ = 1.01
const σ_const = 0.87

f = (t,u) -> μ * u + μ * u
f1 = (t,u) -> μ
(p::typeof(f1))(::Type{Val{:analytic}},t,u0,W) = u0.*exp.((2μ-(σ_const^2)/2)t+σ_const*W)
f2 = (t,u) -> μ * u
σ = (t,u) -> σ_const*u

prob = SplitSDEProblem((f1,f2),σ,1/2,(0.0,1.0))

sol = solve(prob,IIF1M(),dt=1/10)

prob2 = SDEProblem(f,σ,1/2,(0.0,1.0),noise = NoiseWrapper(sol.W))

sol2 = solve(prob2,EM(),dt=1/10)

srand(100)
dts = 1./2.^(7:-1:4) #14->7 good plot

sim  = test_convergence(dts,prob,IIF1M(),numMonte=Int(2e1))
@test abs(sim.𝒪est[:l2]-1) < 0.2 # closer to 1 at this part

srand(200)
sim  = test_convergence(dts,prob,IIF1Mil(),numMonte=Int(2e1))
@test abs(sim.𝒪est[:l2]-1) < 0.3



u0 = rand(2)
A = Strang(2)
B = [σ_const 0
    0 σ_const]

f = function (t,u,du)
  A_mul_B!(du,A,u)
  du .+= 1.01u
end
σ = function (t,u,du)
  A_mul_B!(@view(du[:,1]),B,u)
  A_mul_B!(@view(du[:,2]),B,u)
end

function (p::typeof(f))(::Type{Val{:analytic}},t,u0,W)
 tmp = (A+1.01I-(B^2))*t + B*sum(W)
 expm(tmp)*u0
end

f1 = (t,u,du) -> A
function (p::typeof(f1))(::Type{Val{:analytic}},t,u0,W)
 tmp = (A+1.01I-(B^2))*t + B*sum(W)
 expm(tmp)*u0
end
f2 = (t,u,du) -> du .= μ .* u

prob = SplitSDEProblem((f1,f2),σ,u0,(0.0,1.0),noise_rate_prototype=rand(2,2))

sol = solve(prob,IIF1M(),dt=1/10)

dts = 1./2.^(7:-1:4) #14->7 good plot

srand(250)
sim  = test_convergence(dts,prob,IIF1M(),numMonte=Int(2e1))
@test abs(sim.𝒪est[:l2]-0.5) < 0.2

#=
sim  = test_convergence(dts,prob,IIF1Mil(),numMonte=Int(5e1))
@test abs(sim.𝒪est[:l2]-1) < 0.1
=#
