using DiffEqBase, StochasticDiffEq, Test
explicit = [SRIW1()]
implicit_autodiff   = [SKenCarp(), ImplicitEulerHeun()]
implicit_noautodiff = [SKenCarp(autodiff=false), ImplicitEulerHeun(autodiff=false)]
oopf(u,p,t)=1.01u
oopg(u,p,t)=1.01u
u0 = 1.0+1.0im
tspan = (0.0,1.0)
prob_oop = SDEProblem(oopf,oopg,u0,tspan)
iipf(du,u,p,t) = (du.=1.01u)
iipg(du,u,p,t) = (du.=1.01u)
u0 = ones(2,4) + im*ones(2,4)
prob_iip = SDEProblem(iipf,iipg,u0,tspan)

for alg in vcat(explicit, implicit_noautodiff), prob in [prob_iip, prob_oop]
  sol = solve(prob,alg)
  @test eltype(sol[end]) == ComplexF64
end

for alg in implicit_autodiff, prob in [prob_iip, prob_oop]
  @test_broken sol = solve(prob,alg)
  #eltype(sol[end]) == ComplexF64
end
