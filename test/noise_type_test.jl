using StochasticDiffEq, DiffEqBase, Test, SparseArrays

f(du,u,p,t) = (du.=1.01u)
function g(du,u,p,t)
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

f(du,u,p,t) = (du.=1.01u)
function g(du,u,p,t)
  du[1,1] = 0.3
  du[1,2] = 0.6
  du[1,3] = 0.9
  du[1,4] = 0.12
  du[2,1] = 1.2
  du[2,2] = 0.2
  du[2,3] = 0.3
  du[2,4] = 1.8
end
prob = SDEProblem(f,g,ones(2),(0.0,1.0),noise_rate_prototype=zeros(2,4))

sol = solve(prob,SRA1())
sol = solve(prob,SRA2())
sol = solve(prob,SRA3())
sol = solve(prob,SOSRA())
sol = solve(prob,SOSRA2())
sol = solve(prob,SRA())

@test length(sol.W[1]) == 4

function g(du,u,p,t)
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
