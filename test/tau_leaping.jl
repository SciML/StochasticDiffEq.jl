using StochasticDiffEq, DiffEqJump, DiffEqBase
using Test, LinearAlgebra

function regular_rate(out,u,p,t)
    out[1] = (0.1/1000.0)*u[1]*u[2]
    out[2] = 0.01u[2]
end

const dc = zeros(3, 2)
dc[1,1] = -1
dc[2,1] = 1
dc[2,2] = -1
dc[3,2] = 1

function regular_c(du,u,p,t,counts,mark)
    mul!(du,dc,counts)
end

rj = RegularJump(regular_rate,regular_c,2)
jumps = JumpSet(rj)
prob = DiscreteProblem([999.0,1,0],(0.0,250.0))
jump_prob = JumpProblem(prob,Direct(),rj)
sol = solve(jump_prob,TauLeaping();dt=1.0)
sol = solve(jump_prob,SimpleTauLeaping();dt=1.0)

function rate_oop(u,p,t)
    [(0.1/1000.0)*u[1]*u[2],0.01u[2]]
end

const dc = zeros(3, 2)
dc[1,1] = -1
dc[2,1] = 1
dc[2,2] = -1
dc[3,2] = 1

function regular_c(u,p,t,counts,mark)
    dc*counts
end

rj = RegularJump(rate_oop,regular_c,2)
jumps = JumpSet(rj)
prob = DiscreteProblem([999.0,1,0],(0.0,250.0))
jump_prob = JumpProblem(prob,Direct(),rj)
sol = solve(jump_prob,TauLeaping();dt=1.0)
