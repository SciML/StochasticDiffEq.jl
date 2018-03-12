using DiffEqBase, StochasticDiffEq, DiffEqNoiseProcess, Base.Test, DiffEqDevTools

f(u,p,t) = (1.01) * u
f1(u,p,t) = (1.01)/2 * u
f2(u,p,t) = (1.01)/2 * u
σ(u,p,t) = 0.87u
#(::typeof(f))(::Type{Val{:analytic}},u0,p,t,W) = u0.*exp.(0.63155t+0.87W)

prob = SDEProblem{false}((f1,f2),σ,1/2,(0.0,1.0))

sol = solve(prob,SplitEM(),dt=1/10,save_noise=true)

prob = SDEProblem{false}(f,σ,1/2,(0.0,1.0),noise = NoiseWrapper(sol.W))

sol2 = solve(prob,EM(),dt=1/10)

@test sol[:] ≈ sol2[:]

u0 = rand(4)
prob = SDEProblem{false}((f1,f2),σ,u0,(0.0,1.0))

sol = solve(prob,SplitEM(),dt=1/10,save_noise=true)

prob = SDEProblem{false}(f,σ,u0,(0.0,1.0),noise = NoiseWrapper(sol.W))

sol2 = solve(prob,EM(),dt=1/10)

@test sol[end][:] ≈ sol2[end][:]

################################################################################

### Only first

α = 0.1
β = 0.5
ff1 = (u,p,t) -> β./sqrt.(1+t) - u./(2*(1+t))
ff2 = (u,p,t) -> 0.0
σ2 = (u,p,t) -> α*β./sqrt.(1+t)
prob = SplitSDEProblem(ff1,ff2,σ2,1.,(0.0,1.0))
(::typeof(prob.f))(::Type{Val{:analytic}},u0,p,t,W) = u0./sqrt.(1+t) + β*(t+α*W)./sqrt.(1+t)

sol = solve(prob,EM(),dt=1/10)
sol2 = solve(prob,SKenCarp(),dt=1/10)

dts = 1./2.^(10:-1:2) #14->7 good plot
sim10 = test_convergence(dts,prob,SKenCarp(),numMonte=Int(1e1))
@test abs(sim10.𝒪est[:final]-2) < 0.3

### Only second

α = 0.1
β = 0.5
ff1 = (u,p,t) -> 0.0
ff2 = (u,p,t) -> β./sqrt.(1+t) - u./(2*(1+t))
σ2 = (u,p,t) -> α*β./sqrt.(1+t)
prob = SplitSDEProblem(ff1,ff2,σ2,1.,(0.0,1.0))
(::typeof(prob.f))(::Type{Val{:analytic}},u0,p,t,W) = u0./sqrt.(1+t) + β*(t+α*W)./sqrt.(1+t)

sol = solve(prob,EM(),dt=1/10)
sol2 = solve(prob,SKenCarp(),dt=1/10,seed=1)

dts = 1./2.^(10:-1:2) #14->7 good plot
sim10 = test_convergence(dts,prob,SKenCarp(),numMonte=Int(1e1))
@test abs(sim10.𝒪est[:final]-2) < 0.3

### Both

α = 0.1
β = 0.5
ff1 = (u,p,t) -> β./sqrt.(1+t)
ff2 = (u,p,t) -> - u./(2*(1+t))
σ2 = (u,p,t) -> α*β./sqrt.(1+t)
prob = SplitSDEProblem(ff1,ff2,σ2,1.,(0.0,1.0))
(::typeof(prob.f))(::Type{Val{:analytic}},u0,p,t,W) = u0./sqrt.(1+t) + β*(t+α*W)./sqrt.(1+t)

sol = solve(prob,EM(),dt=1/10)
sol2 = solve(prob,SKenCarp(),dt=1/10)

dts = 1./2.^(10:-1:2) #14->7 good plot
sim10 = test_convergence(dts,prob,SKenCarp(),numMonte=Int(1e1))
@test abs(sim10.𝒪est[:final]-2) < 0.3

################################################################################

### Only first

α = 0.1
β = 0.5
ff1 = (du,u,p,t) -> du .= β./sqrt.(1+t) - u./(2*(1+t))
ff2 = (du,u,p,t) -> du .= 0.0
σ2 = (du,u,p,t) -> du .= α*β./sqrt.(1+t)
prob = SplitSDEProblem(ff1,ff2,σ2,[1.],(0.0,1.0))
(::typeof(prob.f))(::Type{Val{:analytic}},u0,p,t,W) = u0./sqrt.(1+t) + β*(t+α*W)./sqrt.(1+t)

sol = solve(prob,EM(),dt=1/10)
sol2 = solve(prob,SKenCarp(),dt=1/10)

dts = 1./2.^(10:-1:2) #14->7 good plot
sim10 = test_convergence(dts,prob,SKenCarp(),numMonte=Int(1e1))
@test abs(sim10.𝒪est[:final]-2) < 0.3

### Only second

α = 0.1
β = 0.5
ff1 = (du,u,p,t) -> du .= 0.0
ff2 = (du,u,p,t) -> du .= β./sqrt.(1+t) - u./(2*(1+t))
σ2 = (du,u,p,t) -> du .= α*β./sqrt.(1+t)
prob = SplitSDEProblem(ff1,ff2,σ2,[1.],(0.0,1.0))
(::typeof(prob.f))(::Type{Val{:analytic}},u0,p,t,W) = u0./sqrt.(1+t) + β*(t+α*W)./sqrt.(1+t)

sol = solve(prob,EM(),dt=1/10)
sol2 = solve(prob,SKenCarp(),dt=1/10)

dts = 1./2.^(10:-1:2) #14->7 good plot
sim10 = test_convergence(dts,prob,SKenCarp(),numMonte=Int(1e1))
@test abs(sim10.𝒪est[:final]-2) < 0.3

### Both

α = 0.1
β = 0.5
ff1 = (du,u,p,t) -> du .= β./sqrt.(1+t)
ff2 = (du,u,p,t) -> du .= - u./(2*(1+t))
σ2 = (du,u,p,t) -> du .= α*β./sqrt.(1+t)
prob = SplitSDEProblem(ff1,ff2,σ2,[1.],(0.0,1.0))
(::typeof(prob.f))(::Type{Val{:analytic}},u0,p,t,W) = u0./sqrt.(1+t) + β*(t+α*W)./sqrt.(1+t)

sol = solve(prob,EM(),dt=1/10)
sol2 = solve(prob,SKenCarp(),dt=1/10)

dts = 1./2.^(10:-1:2) #14->7 good plot
sim10 = test_convergence(dts,prob,SKenCarp(),numMonte=Int(1e1))
@test abs(sim10.𝒪est[:final]-2) < 0.3
