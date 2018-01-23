using DiffEqBase, StochasticDiffEq, DiffEqDevTools, Base.Test

u0 = ones(2)
A = [-3/2 1/20
      1/20 -3/2]
B = [1/5 1/100
    1/100 1/5]

function f(du,u,p,t)
  A_mul_B!(du,A,u)
end
function σ(du,u,p,t)
  A_mul_B!(@view(du[:,1]),B,u)
  A_mul_B!(@view(du[:,2]),B,u)
end

function (::typeof(f))(::Type{Val{:analytic}},u0,p,t,W)
  tmp = (A-(B^2))*t + B*W[1] + B*W[2]
  expm(tmp)*u0
end

prob2 = SDEProblem(f,σ,u0,(0.0,1.0),noise_rate_prototype=rand(2,2))

sol2 = solve(prob2,EM(),dt=1/100)
using Plots; plot(sol2,plot_analytic=true)

dts = 1./2.^(14:-1:7) #14->7 good plot

sim  = test_convergence(dts,prob2,EM(),numMonte=Int(5e1))
@test abs(sim.𝒪est[:l2]-0.5) < 0.1

using Plots; plot(sim)


using SpecialMatrices
u0 = rand(2)
A = Strang(2)
B = [1/5 1/100
    1/100 1/5]

function f(du,u,p,t)
  A_mul_B!(du,A,u)
  du .+= 1.01u
end
function σ(du,u,p,t)
  A_mul_B!(@view(du[:,1]),B,u)
  A_mul_B!(@view(du[:,2]),B,u)
end

function (::typeof(f))(::Type{Val{:analytic}},u0,p,t,W)
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
