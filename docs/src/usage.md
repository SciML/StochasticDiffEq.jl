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
sol = solve(prob, EM(split=true))

# RKMilCommute with Stratonovich interpretation
sol = solve(prob, RKMilCommute(interpretation=:Stratonovich))

# Implicit methods with solver options
sol = solve(prob, SKenCarp(linsolve=KrylovJL_GMRES()))
```

## Tolerances and Adaptive Stepping

Set absolute and relative tolerances:

```julia
sol = solve(prob, SOSRI(), abstol=1e-6, reltol=1e-3)
```

For fixed time stepping:
```julia  
sol = solve(prob, EM(), dt=0.01, adaptive=false)
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

### Using the Same Noise Realization Multiple Times

For applications like computing snapshot attractors of random dynamical systems, you may need to solve the same SDE with different initial conditions but the same noise realization. This requires careful handling of the random number generator and noise process.

#### Basic Approach

```julia
using StochasticDiffEq, Random, DiffEqNoiseProcess

# Define your SDE
function f!(du, u, p, t)
    du[1] = -u[1] + p[1] * u[2]
    du[2] = -u[2] - p[1] * u[1]
end

function g!(du, u, p, t)
    du[1] = 0.1
    du[2] = 0.1
end

# Create a specific RNG with a non-zero seed
rng = Random.Xoshiro(1234)

# Create a noise process that won't reseed
noise = WienerProcess(0.0, 0.0, rng=rng, reseed=false)

p = [2.0]
tspan = (0.0, 1.0)

# Solve with first initial condition
u0_1 = [1.0, 0.0]
prob1 = SDEProblem(f!, g!, u0_1, tspan, p, noise=noise)
sol1 = solve(prob1, EM(), dt=0.01)

# Reset the noise process to reuse the same realization
# Create a new noise process with the same seed
rng_reset = Random.Xoshiro(1234)  # Same seed as before
noise = WienerProcess(0.0, 0.0, rng=rng_reset, reseed=false)

# Solve with second initial condition using the same noise
u0_2 = [0.0, 1.0]
prob2 = SDEProblem(f!, g!, u0_2, tspan, p, noise=noise)
sol2 = solve(prob2, EM(), dt=0.01)

# Verify both solutions used the same noise
@assert sol1.W.W == sol2.W.W  # Same noise realization
```

#### Solving Across Time Windows

To solve an SDE across different time windows with the same noise continuation:

```julia
using StochasticDiffEq, Random, DiffEqNoiseProcess

# Define SDE (same as above)
function f!(du, u, p, t)
    du[1] = -0.1 * u[1]
    du[2] = -0.1 * u[2]
end

function g!(du, u, p, t)
    du[1] = 0.2
    du[2] = 0.2
end

# Setup with persistent RNG
rng = Random.Xoshiro(5678)
noise = WienerProcess(0.0, 0.0, rng=rng, reseed=false)

u0 = [1.0, 1.0]
dt = 0.01

# First time window: t ∈ [0, 1]
prob1 = SDEProblem(f!, g!, u0, (0.0, 1.0), noise=noise)
sol1 = solve(prob1, EM(), dt=dt)

# Second time window: t ∈ [1, 2] 
# Continue from the end state of the first solution
# The noise process automatically continues from where it left off
u0_2 = sol1.u[end]  # Final state from first window
prob2 = SDEProblem(f!, g!, u0_2, (1.0, 2.0), noise=noise)
sol2 = solve(prob2, EM(), dt=dt)

# Combine solutions if needed
combined_t = [sol1.t[1:end-1]; sol2.t]  # Avoid duplicate at t=1
combined_u = [sol1.u[1:end-1]; sol2.u]
```

#### Important Notes

1. **Use non-zero seed**: Never use `seed=0` as this can cause issues with some RNGs
2. **Set `reseed=false`**: This prevents the noise process from regenerating 
3. **Consistent timestep**: Use the same `dt` across different solves for identical noise
4. **Reset when needed**: Create a new noise process with the same seed to restart the same noise realization
5. **Compatible solvers**: This approach works with fixed-timestep methods like `EM()`

This technique is particularly valuable for:
- Computing snapshot attractors in random dynamical systems
- Studying sensitivity to initial conditions with fixed noise
- Reproducible ensemble simulations with shared noise components

## Itô vs Stratonovich

Specify interpretation when creating problems or choosing solvers:

```julia
# Itô interpretation (default)
prob = SDEProblem(f!, g!, u0, tspan, interpretation=:Ito)

# Stratonovich interpretation  
prob = SDEProblem(f!, g!, u0, tspan, interpretation=:Stratonovich)

# Or at solver level
sol = solve(prob, RKMil(interpretation=:Stratonovich))
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
condition(u,t,integrator) = u[1] - 0.5
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition, affect!)

sol = solve(prob, SOSRI(), callback=cb)

# Ensemble simulations
monte_prob = EnsembleProblem(prob)
sim = solve(monte_prob, SOSRI(), EnsembleThreads(), trajectories=1000)
```