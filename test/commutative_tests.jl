@everywhere using StochasticDiffEq, DiffEqProblemLibrary, DiffEqDevTools, Base.Test
srand(100)
dts = 1./2.^(10:-1:2) #14->7 good plot

using SpecialMatrices
const σ_const = 0.87

u0 = rand(2)
A = Strang(2)
B = [σ_const 0
    0 σ_const]

function f(t,u,du)
  A_mul_B!(du,A,u)
end
function (p::typeof(f))(::Type{Val{:analytic}},t,u0,W)
 tmp = (A-(B^2))*t + B*sum(W)
 expm(tmp)*u0
end
function σ(t,u,du)
  du[1,1] = σ_const*u[1]
  du[1,2] = σ_const*u[1]
  du[2,1] = σ_const*u[2]
  du[2,2] = σ_const*u[2]
end
g = function (t,u,du)
  du .= σ_const.*u
end

prob = SDEProblem(f,σ,u0,(0.0,1.0),noise_rate_prototype=rand(2,2))

sol = solve(prob,RKMilCommute(),dt=1/2^(8))
sol = solve(prob,EM(),dt=1/2^(10))

dts = 1./2.^(10:-1:3) #14->7 good plot
sim2 = test_convergence(dts,prob,RKMilCommute(),numMonte=Int(1e2))












sim2 = test_convergence(dts,prob,EM(),numMonte=Int(1e2))











@test abs(sim2.𝒪est[:l∞]-1) < 0.2

dW = rand(2); dt = 1/2; uprev = u0; t = 0; du1 = similar(u0); sqdt = sqrt(dt)
I = zeros(length(dW),length(dW)); L = zeros(length(dW),length(dW)); u = similar(uprev)
Dg = zeros(length(dW),length(dW)); mil_correction = zeros(length(dW))
tmp = similar(L)
mil_correction .= 0.0
for i=1:length(dW),j=1:length(dW)
    I[j,i] = 0.5*dW[i]*dW[j]
    j == i && (I[i,i] -= 0.5*dt) # Ito correction
end

I[1,1] == (dW[1].^2 - dt)/2

f(t,uprev,du1)


σ(t,uprev,L)

for j = 1:length(uprev)
  Kj = uprev .+ dt.*du1 + sqdt*L[j,:]
  σ(t,Kj,tmp)
  Dgj = (tmp - L)/sqdt
  mil_correction .+= Dgj*I[j,:]
end
u .= uprev .+ dt.*du1 + L*dW# .+ mil_correction




σ = function (t,u,du)
  A_mul_B!(@view(du[:,1]),B,u)
  A_mul_B!(@view(du[:,2]),B,u)
end
function (p::typeof(f))(::Type{Val{:analytic}},t,u0,W)
 tmp = (A+1.01I-(B^2))*t + B*sum(W)
 expm(tmp)*u0
end
