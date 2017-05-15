using StaticArrays
using DiffEqBase, StochasticDiffEq

#=
u0 = zeros(MVector{2,Float64}, 2) + 1
u0[1] = ones(MVector{2,Float64}) + 1
f = (t,u,du) -> du .= u
ode = SDEProblem(f,f, u0, (0.,1.))
sol = solve(ode, EM(), dt=1.e-2)

u0 = zeros(SVector{2,Float64}, 2) + 1
u0[1] = ones(SVector{2,Float64}) + 1
ode = SDEProblem(f, f, u0, (0.,1.))
sol = solve(ode, EM(), dt=1.e-2)
sol = solve(ode, SRIW1(), dt=1.e-2)
=#

u0 = zero(MVector{2,Float64}) + 1
ode = SDEProblem(f, f, u0, (0.,1.))
sol = solve(ode, EM(), dt=1.e-2)
sol = solve(ode, SRIW1())

u0 = zero(SVector{2,Float64}) + 1
f = (t,u) -> u
prob = SDEProblem(f, f, u0, (0.,1.))
sol = solve(ode, EM(), dt=1.e-2)
sol = solve(ode, SRIW1())
