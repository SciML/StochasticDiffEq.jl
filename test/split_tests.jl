using DiffEqBase, StochasticDiffEq, DiffEqNoiseProcess, Base.Test

f = (t,u) -> (1.01) * u
f1 = (t,u) -> (1.01)/2 * u
f2 = (t,u) -> (1.01)/2 * u
σ = (t,u) -> 0.87u
#(p::typeof(f))(::Type{Val{:analytic}},t,u0,W) = u0.*exp.(0.63155t+0.87W)

prob = SplitSDEProblem((f1,f2),σ,1/2,(0.0,1.0))

sol = solve(prob,SplitEM(),dt=1/10)

prob = SDEProblem(f,σ,1/2,(0.0,1.0),noise = NoiseWrapper(sol.W))

sol2 = solve(prob,EM(),dt=1/10)

@test sol[:] ≈ sol2[:]

u0 = rand(4)
prob = SplitSDEProblem((f1,f2),σ,u0,(0.0,1.0))

sol = solve(prob,SplitEM(),dt=1/10)

prob = SDEProblem(f,σ,u0,(0.0,1.0),noise = NoiseWrapper(sol.W))

sol2 = solve(prob,EM(),dt=1/10)

@test sol[end][:] ≈ sol2[end][:]
