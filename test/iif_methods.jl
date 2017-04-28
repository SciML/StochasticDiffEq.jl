using DiffEqBase, StochasticDiffEq, DiffEqNoiseProcess, Base.Test, DiffEqDevTools

using SpecialMatrices

const μ = 1.01
const σ_const = 0.87
f = (t,u) -> A * u + μ * u
f1 = (t,u) -> A
(p::typeof(f1))(::Type{Val{:analytic}},t,u0,W) = u0.*exp.(0.63155t+σ_const*W)
f2 = (t,u) -> μ * u
σ = (t,u) -> σ_const*u
#(p::typeof(f))(::Type{Val{:analytic}},t,u0,W) = u0.*exp.(0.63155t+σ_constW)

prob = SplitSDEProblem((f1,f2),σ,1/2,(0.0,1.0))

sol = solve(prob,IIF1M(),dt=1/10000)

prob2 = SDEProblem(f,σ,1/2,(0.0,1.0),noise = NoiseWrapper(sol.W))

sol2 = solve(prob2,EM(),dt=1/10000)

srand(100)
dts = 1./2.^(10:-1:2) #14->7 good plot

sim  = test_convergence(dts,prob,IIF1M(),numMonte=Int(1e3))
@test abs(sim.𝒪est[:l2]-0.5) < 0.1

sim  = test_convergence(dts,prob,IIF1Mil(),numMonte=Int(1e3))
@test abs(sim.𝒪est[:l2]-1) < 0.1



using SpecialMatrices
u0 = rand(2)
A = Strang(2)
B = [1/5 1/100
    1/100 1/5]

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

prob2 = SDEProblem(f,σ,u0,(0.0,1.0),noise_rate_prototype=rand(2,2))

sol2 = solve(prob2,EM(),dt=1/100)
using Plots; plot(sol2,plot_analytic=true)

dts = 1./2.^(17:-1:10) #14->7 good plot

sim  = test_convergence(dts,prob2,EM(),numMonte=Int(5e1))
@test abs(sim.𝒪est[:l2]-0.5) < 0.1

using Plots; plot(sim)


f1 = (t,u,du) -> A
function (p::typeof(f1))(::Type{Val{:analytic}},t,u0,W)
 tmp = (A+1.01I-(B^2))*t + B*sum(W)
 expm(tmp)*u0
end
f2 = (t,u,du) -> du .= μ .* u

prob = SplitSDEProblem((f1,f2),σ,u0,(0.0,1.0))

sol = solve(prob,IIF1M(),dt=1/10000)

sim  = test_convergence(dts,prob,IIF1M(),numMonte=Int(1e3))
@test abs(sim.𝒪est[:l2]-0.5) < 0.1

sim  = test_convergence(dts,prob,IIF1Mil(),numMonte=Int(1e3))
@test abs(sim.𝒪est[:l2]-1) < 0.1
