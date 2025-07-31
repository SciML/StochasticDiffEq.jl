# Usage Guide

This page provides guidance on using StochasticDiffEq.jl effectively.

## Basic Usage

### Problem Definition

StochasticDiffEq.jl uses the standard DifferentialEquations.jl problem interface:

```julia
using StochasticDiffEq

# For scalar problems
function f(u, p, t)  # drift
    return μ * u
end

function g(u, p, t)  # diffusion  
    return σ * u
end

# For in-place systems
function f!(du, u, p, t)
    du[1] = μ * u[1]
    du[2] = -ν * u[2]
end

function g!(du, u, p, t)
    du[1] = σ₁ * u[1]
    du[2] = σ₂ * u[2]
end

# Create problem
prob = SDEProblem(f, g, u0, tspan)
```

### Solver Selection

Choose solvers based on your problem characteristics:

```julia
# Default - good for most problems
sol = solve(prob)

# Specify solver explicitly
sol = solve(prob, SOSRI())          # Recommended for diagonal noise
sol = solve(prob, SOSRA())          # Optimal for additive noise  
sol = solve(prob, SKenCarp())       # For stiff problems
sol = solve(prob, EM())             # For maximum efficiency
```

### Algorithm Parameters

Most solvers accept parameters for customization:

```julia
# Euler-Maruyama with step splitting
sol = solve(prob, EM(split = true))

# RKMilCommute with Stratonovich interpretation
sol = solve(prob, RKMilCommute(interpretation = :Stratonovich))

# Implicit methods with solver options
sol = solve(prob, SKenCarp(linsolve = KrylovJL_GMRES()))
```

## Tolerances and Adaptive Stepping

Set absolute and relative tolerances:

```julia
sol = solve(prob, SOSRI(), abstol = 1e-6, reltol = 1e-3)
```

For fixed time stepping:

```julia
sol = solve(prob, EM(), dt = 0.01, adaptive = false)
```

## Noise Types

### Diagonal Noise

Most common case - each component has independent noise:

```julia
function g!(du, u, p, t)
    du[1] = σ₁ * u[1]
    du[2] = σ₂ * u[2]
end
```

### Scalar Noise

Single noise source affects all components:

```julia
function g!(du, u, p, t)
    du[1] = σ * u[1]
    du[2] = σ * u[2]
end
```

### Non-diagonal Noise

Multiple noise sources with cross-terms:

```julia
function g!(du, u, p, t)
    du[1] = σ₁₁ * u[1] + σ₁₂ * u[2]
    du[2] = σ₂₁ * u[1] + σ₂₂ * u[2]
end
```

### Additive Noise

Noise independent of solution:

```julia
function g!(du, u, p, t)
    du[1] = σ₁
    du[2] = σ₂
end
```

## Itô vs Stratonovich

Specify interpretation when creating problems or choosing solvers:

```julia
# Itô interpretation (default)
prob = SDEProblem(f!, g!, u0, tspan, interpretation = :Ito)

# Stratonovich interpretation  
prob = SDEProblem(f!, g!, u0, tspan, interpretation = :Stratonovich)

# Or at solver level
sol = solve(prob, RKMil(interpretation = :Stratonovich))
```

## Performance Tips

 1. **Use appropriate solvers**: Match solver to problem type
 2. **In-place functions**: Use `f!(du,u,p,t)` for better performance
 3. **Tolerances**: Don't make tolerances unnecessarily strict
 4. **Static arrays**: Use `StaticArrays.jl` for small systems
 5. **GPU**: Use `CuArrays.jl` for large problems

## Common Pitfalls

 1. **Wrong noise type**: Ensure solver supports your noise structure
 2. **Stiffness**: Use appropriate stiff solvers for stiff problems
 3. **Commuting noise**: Use specialized solvers for better efficiency
 4. **High dimensions**: Consider weak convergence methods for Monte Carlo

## Integration with DifferentialEquations.jl

StochasticDiffEq.jl integrates with the broader ecosystem:

```julia
using DifferentialEquations

# Callbacks
condition(u, t, integrator) = u[1] - 0.5
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition, affect!)

sol = solve(prob, SOSRI(), callback = cb)

# Ensemble simulations
monte_prob = EnsembleProblem(prob)
sim = solve(monte_prob, SOSRI(), EnsembleThreads(), trajectories = 1000)
```
