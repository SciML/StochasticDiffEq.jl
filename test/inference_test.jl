using StochasticDiffEq, ForwardDiff
using Test

u0 = 1 / 2
f(u, p, t) = u
g(u, p, t) = u
dt = 1 // 2^(4)
tspan = (0.0, 1.0)
prob = SDEProblem{false}(f, g, u0, (0.0, 1.0))
sol = @inferred solve(prob, EM(), dt = dt)

# issue #351
@inferred StochasticDiffEq.sde_interpolation(ForwardDiff.Dual(0.5), sol.interp, 1, Val{0},
                                             prob.p, :left)
