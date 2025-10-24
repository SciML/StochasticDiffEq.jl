using StochasticDiffEq, JumpProcesses, DiffEqBase, Statistics, OrdinaryDiffEq
using Test, LinearAlgebra

function regular_rate(out, u, p, t)
    out[1] = (0.1/1000.0)*u[1]*u[2]
    out[2] = 0.01u[2]
end

const dc = zeros(3, 2)
dc[1, 1] = -1
dc[2, 1] = 1
dc[2, 2] = -1
dc[3, 2] = 1

function regular_c(du, u, p, t, counts, mark)
    mul!(du, dc, counts)
end

rj = RegularJump(regular_rate, regular_c, 2)
jumps = JumpSet(rj)
iip_prob = DiscreteProblem([999.0, 1, 0], (0.0, 250.0))
jump_iipprob = JumpProblem(iip_prob, Direct(), rj)
@time sol = solve(jump_iipprob, TauLeaping())
@time sol = solve(jump_iipprob, SimpleTauLeaping(); dt = 1.0)
@time sol = solve(jump_iipprob, TauLeaping(); dt = 1.0, adaptive = false)
@time sol = solve(jump_iipprob, CaoTauLeaping(); dt = 1.0)
@time sol = solve(jump_iipprob, CaoTauLeaping())

N = 40_000
sol1 = solve(EnsembleProblem(jump_iipprob), SimpleTauLeaping(); dt = 1.0, trajectories = N)
sol2 = solve(EnsembleProblem(jump_iipprob), TauLeaping(); dt = 1.0,
    adaptive = false, save_everystep = false, trajectories = N)

mean1 = mean([sol1.u[i][end, end] for i in 1:N])
mean2 = mean([sol2.u[i][end, end] for i in 1:N])
@test mean1 ≈ mean2 rtol=1e-2

f(du, u, p, t) = (du .= 0)
g(du, u, p, t) = (du .= 0)
iip_sdeprob = SDEProblem(f, g, [999.0, 1, 0], (0.0, 250.0))
jumpdiff_iipprob = JumpProblem(iip_sdeprob, Direct(), rj)
@time sol = solve(jumpdiff_iipprob, EM(); dt = 1.0)
@time sol = solve(jumpdiff_iipprob, ImplicitEM(); dt = 1.0, adaptive = false)

sol = solve(EnsembleProblem(jumpdiff_iipprob), EM(); dt = 1.0, trajectories = N)
meanX = mean([sol.u[i][end, end] for i in 1:N])
@test mean1 ≈ meanX rtol=1e-2

sol = solve(EnsembleProblem(jumpdiff_iipprob), ImplicitEM(); dt = 1.0, trajectories = N)
meanX = mean([sol.u[i][end, end] for i in 1:N])
@test mean1 ≈ meanX rtol=1e-2

iip_prob = DiscreteProblem([999, 1, 0], (0.0, 250.0))
jump_iipprob = JumpProblem(iip_prob, Direct(), rj)
sol = solve(jump_iipprob, TauLeaping())

function rate_oop(u, p, t)
    [(0.1/1000.0)*u[1]*u[2], 0.01u[2]]
end

function regular_c(u, p, t, counts, mark)
    dc*counts
end

rj = RegularJump(rate_oop, regular_c, 2)
jumps = JumpSet(rj)
prob = DiscreteProblem([999.0, 1, 0], (0.0, 250.0))
jump_prob = JumpProblem(prob, Direct(), rj)
sol = solve(jump_prob, TauLeaping(), reltol = 5e-2)

sol2 = solve(EnsembleProblem(jump_prob), TauLeaping(); dt = 1.0,
    adaptive = false, save_everystep = false, trajectories = N)
mean2 = mean([sol2.u[i][end, end] for i in 1:N])
@test mean1 ≈ mean2 rtol=1e-2

foop(u, p, t) = [0.0, 0.0, 0.0]
goop(u, p, t) = [0.0, 0.0, 0.0]
oop_sdeprob = SDEProblem(foop, goop, [999.0, 1, 0], (0.0, 250.0))
jumpdiff_prob = JumpProblem(oop_sdeprob, Direct(), rj)
@time sol = solve(jumpdiff_prob, EM(); dt = 1.0)
@time sol = solve(jumpdiff_prob, ImplicitEM(); dt = 1.0)

sol = solve(EnsembleProblem(jumpdiff_prob), EM(); dt = 1.0, trajectories = 10_000)
meanX = mean([sol.u[i][end, end] for i in 1:10_000])
@test mean1 ≈ meanX rtol=1e-2

sol = solve(EnsembleProblem(jumpdiff_prob), ImplicitEM(); dt = 1.0, trajectories = 1_000)
meanX = mean([sol.u[i][end, end] for i in 1:1_000])
@test mean1 ≈ meanX rtol=1e-1
