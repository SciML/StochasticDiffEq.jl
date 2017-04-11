using StochasticDiffEq, DiffEqBase, Base.Test

f = (t,u,du) -> du.=1.01u
g = function (t,u,du)
  du[1,1] = 0.3u[1]
  du[1,2] = 0.6u[1]
  du[1,3] = 0.9u[1]
  du[1,4] = 0.12u[2]
  du[2,1] = 1.2u[1]
  du[2,2] = 0.2u[2]
  du[2,3] = 0.3u[2]
  du[2,4] = 1.8u[2]
end
prob = SDEProblem(f,g,ones(2),(0.0,1.0),noise_rate_prototype=zeros(2,4))

sol = solve(prob,EM(),dt=1/1000)

@test length(sol.W[1]) == 4

g = function (t,u,du)
  @test typeof(du) <: SparseMatrixCSC
  du[1,1] = 0.3u[1]
  du[1,2] = 0.6u[1]
  du[1,3] = 0.9u[1]
  du[1,4] = 0.12u[2]
  du[2,1] = 1.2u[1]
  du[2,2] = 0.2u[2]
  du[2,3] = 0.3u[2]
  du[2,4] = 1.8u[2]
end
prob = SDEProblem(f,g,ones(2),(0.0,1.0),noise_rate_prototype=sprand(2,4,1.0))

sol = solve(prob,EM(),dt=1/1000)

@test length(sol.W[1]) == 4
