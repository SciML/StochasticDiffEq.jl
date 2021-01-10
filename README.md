# StochasticDiffEq.jl 

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://github.com/SciML/StochasticDiffEq.jl/workflows/CI/badge.svg)](https://github.com/SciML/StochasticDiffEq.jl/actions?query=workflow%3ACI)
[![GitlabCI](https://gitlab.com/JuliaGPU/StochasticDiffEq.jl/badges/master/pipeline.svg)](https://gitlab.com/JuliaGPU/StochasticDiffEq.jl/pipelines)

StochasticDiffEq.jl is a component package in the DifferentialEquations ecosystem. It holds the
stochastic differential equations solvers and utilities. While completely independent
and usable on its own, users interested in using this
functionality should check out [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

## API

StochasticDiffEq.jl is part of the JuliaDiffEq common interface, but can be used independently of DifferentialEquations.jl. The only requirement is that the user passes an StochasticDiffEq.jl algorithm to `solve`. For example, we can solve the [SDE tutorial from the docs](https://diffeq.sciml.ai/stable/tutorials/sde_example/) using the `SRIW1()` algorithm:

```julia
using StochasticDiffEq
α=1
β=1
u₀=1/2
f(u,p,t) = α*u
g(u,p,t) = β*u
dt = 1//2^(4)
tspan = (0.0,1.0)
prob = SDEProblem(f,g,u₀,(0.0,1.0))
sol =solve(prob,SRIW1())
```

The options for `solve` are defined in the [common solver options page](https://diffeq.sciml.ai/stable/basics/common_solver_opts/) and are thoroughly explained in [the ODE tutorial](https://diffeq.sciml.ai/stable/tutorials/ode_example/).

That example uses the out-of-place syntax `f(u,p,t)`, while the inplace syntax (more efficient for systems of equations) is shown in the Lorenz example:

```julia
function lorenz(du,u,p,t)
 du[1] = 10.0(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end

function σ_lorenz(du,u,p,t)
 du[1] = 3.0
 du[2] = 3.0
 du[3] = 3.0
end

prob_sde_lorenz = SDEProblem(lorenz,σ_lorenz,[1.0,0.0,0.0],(0.0,10.0))
sol = solve(prob_sde_lorenz)
plot(sol,vars=(1,2,3))
```

The problems default to diagonal noise. Non-diagonal noise can be added by setting
the `noise_prototype`:

```julia
f = (du,u,p,t) -> du.=1.01u
g = function (du,u,p,t)
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
```

Colored noise can be set using [an `AbstractNoiseProcess`](https://diffeq.sciml.ai/stable/features/noise_process/). For example, we can set the underlying noise process to a `GeometricBrownian` via:

```julia
μ = 1.0
σ = 2.0
W = GeometricBrownianMotionProcess(μ,σ,0.0,1.0,1.0)
# ...
# Define f,g,u0,tspan for a SDEProblem
# ...
prob = SDEProblem(f,g,u0,tspan,noise=W)
```

StochasticDiffEq.jl also handles solving random ordinary differential equations. This is shown [in the RODE tutorial](https://diffeq.sciml.ai/stable/tutorials/rode_example/).

```julia
using StochasticDiffEq
function f(u,p,t,W)
  2u*sin(W)
end
u0 = 1.00
tspan = (0.0,5.0)
prob = RODEProblem(f,u0,tspan)
sol = solve(prob,RandomEM(),dt=1/100)
```

## Available Solvers

For the list of available solvers, please refer to the [DifferentialEquations.jl SDE Solvers page](https://diffeq.sciml.ai/stable/solvers/sde_solve/) and the [RODE Solvers page](https://diffeq.sciml.ai/stable/solvers/rode_solve/).
