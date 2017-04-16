using DiffEqBase, StochasticDiffEq


type MyNoise{T}
  count::Int
  rands::T
end
function (p::MyNoise)(integrator)
  p.count += 1
  p.rands[p.count]
end
noise = MyNoise(0,randn(100))
my_noise = NoiseProcess{:White,false,typeof(noise)}(noise)


f1 = (t,u) -> 0.
g1 = (t,u) -> 1.
dt = 1//2^(4)
prob1 = SDEProblem(f1,g1,0.,(0.0,1.0),noise=my_noise)
sol1 = solve(prob1,EM(),dt=dt)
noise.count == length(sol1)-1

noise.count = 0
sol1 = solve(prob1,SRIW1(),dt=dt,adaptive=false)
noise.count == 2length(sol1)-2
