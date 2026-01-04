using StochasticDiffEq, DiffEqNoiseProcess, ForwardDiff
using Test

u0 = 1 / 2
f(u, p, t) = u
g(u, p, t) = u
dt = 1 // 2^(4)
tspan = (0.0, 1.0)
prob = SDEProblem{false}(f, g, u0, (0.0, 1.0))
sol = @inferred solve(prob, EM(), dt = dt)

# issue #351
@inferred StochasticDiffEq.sde_interpolation(
    ForwardDiff.Dual(0.5), sol.interp, 1, Val{0}, prob.p, :left
)

# https://github.com/SciML/DiffEqNoiseProcess.jl/issues/213
f(u, p, t) = -u
g(u, p, t) = 1.0

w = WienerProcess(0.0, 0.0)
prob = SDEProblem(f, g, 0.0, (0.0, 1.0), noise = w)
@test_throws "Higher order solver requires extra Brownian process Z. Thus `WienerProcess(t, W0)` is insufficient, you must use `WienerProcess(t, W0, Z0)` where `Z` is another Brownian process" solve(
    prob, SRIW1()
)

w = WienerProcess(0.0, 0.0, 0.0)
prob = SDEProblem(f, g, 0.0, (0.0, 1.0), noise = w)
@test SciMLBase.successful_retcode(solve(prob, SRIW1()))
