using StaticArrays, StochasticDiffEq
using Test
f(du,u,p,t) = (du .= u)

u0 = zeros(MVector{2,Float64}, 2) .+ 1
u0[1] = ones(MVector{2,Float64}) .+ 1
ode = SDEProblem(f,f, u0, (0.,1.))
@test_broken sol = solve(ode, EM(), dt=1.e-2)
@test_broken sol = solve(ode, SRIW1())

u0 = zeros(SVector{2,Float64}, 2) .+ 1
u0[1] = ones(SVector{2,Float64}) .+ 1
ode = SDEProblem(f, f, u0, (0.,1.))
@test_broken sol = solve(ode, EM(), dt=1.e-2)
@test_broken sol = solve(ode, SRIW1(), dt=1.e-2)

u0 = ones(MVector{2,Float64})
ode = SDEProblem(f, f, u0, (0.,1.))
sol = solve(ode, EM(), dt=1.e-2)
sol = solve(ode, SRIW1())
sol = solve(ode, DRI1(), dt=1.e-2)

u0 = ones(SVector{2,Float64})
f(u,p,t) = u
prob = SDEProblem(f, f, u0, (0.,1.))
@test_broken sol = solve(prob, EM(), dt=1.e-2)
@test_broken sol = solve(prob, SRIW1())
@test_broken sol = solve(prob, DRI1(), dt=1.e-2)
