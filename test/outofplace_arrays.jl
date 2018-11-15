using StochasticDiffEq, DiffEqNoiseProcess

@test_broken begin
  # Need to internally broadcast!
  oup = OrnsteinUhlenbeckProcess!(1.0, 0.0, 1.0, 0.0, [0.0], [0.0])
  noise = NoiseWrapper(oup)
  prob = SDEProblem((u,p,t) -> -u, (u,p,t) -> [1.0], [0.0], (0.0, 1.0), noise=noise)
  sol = solve(prob,SRIW1())
  oup = OrnsteinUhlenbeckProcess!(1.0, 0.0, 1.0, 0.0, [0.0], [0.0])
  noise = NoiseWrapper(oup)
  prob = SDEProblem((u,p,t) -> -u, (u,p,t) -> [1.0], [0.0], (0.0, 1.0), noise=noise)
  sol = solve(prob,SOSRI())
  oup = OrnsteinUhlenbeckProcess!(1.0, 0.0, 1.0, 0.0, [0.0], [0.0])
  noise = NoiseWrapper(oup)
  prob = SDEProblem((u,p,t) -> -u, (u,p,t) -> [1.0], [0.0], (0.0, 1.0), noise=noise)
  @test_broken sol = solve(prob,SRA1())
  oup = OrnsteinUhlenbeckProcess!(1.0, 0.0, 1.0, 0.0, [0.0], [0.0])
  noise = NoiseWrapper(oup)
  prob = SDEProblem((u,p,t) -> -u, (u,p,t) -> [1.0], [0.0], (0.0, 1.0), noise=noise)
  @test_broken sol = solve(prob,SRA2())
  oup = OrnsteinUhlenbeckProcess!(1.0, 0.0, 1.0, 0.0, [0.0], [0.0])
  noise = NoiseWrapper(oup)
  prob = SDEProblem((u,p,t) -> -u, (u,p,t) -> [1.0], [0.0], (0.0, 1.0), noise=noise)
  @test_broken sol = solve(prob,SRA3())

  println("Scalar g")
  A = [-1.0 0.0; 0.0 -0.5]
  u0 = [1.0, 1.0]; tspan = (0.0,1.0)
  _f = (u,p,t) -> t*(A*u)
  _g = (u,p,t) -> 1.0
  prob = SDEProblem(SDEFunction(_f, _g), _g, u0, tspan)
  integrator = init(prob, SKenCarp(); adaptive=false, dt=0.01)
  step!(integrator)

  println("Vector g")
  _g = (u,p,t) -> [1.0, 1.0]
  prob = SDEProblem(SDEFunction(_f, _g), _g, u0, tspan)
  println("Implicit EM")
  integrator = init(prob, ImplicitEM(); adaptive=false, dt=0.01)
  step!(integrator)
  println("SKenCarp")
  integrator = init(prob, SKenCarp(); adaptive=false, dt=0.01)
  step!(integrator)
  solve(prob, SOSRI(); adaptive=false, dt=0.01)
  solve(prob, SOSRA(); adaptive=false, dt=0.01)
  solve(prob, EM(); adaptive=false, dt=0.01)
  solve(prob, RKMil(); adaptive=false, dt=0.01)
  solve(prob, SOSRI2(); adaptive=false, dt=0.01)
  solve(prob, SOSRA2(); adaptive=false, dt=0.01)
end
