function _get_alias_noise_from_kwargs(; alias_noise = nothing, alias = nothing, kwargs...)
    if alias_noise !== nothing
        return alias_noise
    elseif alias !== nothing && hasproperty(alias, :alias_noise) && alias.alias_noise !== nothing
        return alias.alias_noise
    else
        return true
    end
end

function DiffEqBase.__solve(
        prob::DiffEqBase.AbstractRODEProblem,
        alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm};
        kwargs...
    )
    integrator = DiffEqBase.__init(prob, alg; kwargs...)
    solve!(integrator)
    if prob isa DiffEqBase.AbstractRODEProblem &&
            typeof(prob.noise) == typeof(integrator.sol.W) &&
            _get_alias_noise_from_kwargs(; kwargs...)
        copy!(prob.noise, integrator.sol.W)
    end
    return integrator.sol
end

# More specific method for JumpProblem to win over JumpProcesses.jl's ambiguity fix dispatch
function DiffEqBase.__solve(
        prob::JumpProblem,
        alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm};
        merge_callbacks = true, kwargs...
    )
    kwargs = DiffEqBase.merge_problem_kwargs(prob; merge_callbacks, kwargs...)
    integrator = DiffEqBase.__init(prob, alg; kwargs...)
    solve!(integrator)
    if concrete_prob(prob) isa DiffEqBase.AbstractRODEProblem &&
            typeof(concrete_prob(prob).noise) == typeof(integrator.sol.W) &&
            _get_alias_noise_from_kwargs(; kwargs...)
        copy!(concrete_prob(prob).noise, integrator.sol.W)
    end
    return integrator.sol
end

# Make it easy to grab the RODEProblem/SDEProblem/DiscreteProblem from the keyword arguments
concrete_prob(prob) = prob
concrete_prob(prob::JumpProblem) = prob.prob

"""
    _resolve_rng(rng, seed, prob) -> (rng, seed, rng_provided)

Resolve the RNG and seed for an SDE/RODE integration from the user-provided
`rng` and `seed` kwargs plus the problem's stored seed.

## Priority

1. **`rng` provided**: Use it directly. `rng_provided = true`. If the RNG is a
   `TaskLocalRNG`, it is converted to a concrete `Xoshiro` seeded from a draw
   from the task-local stream, and the derived seed is returned as `seed`. This
   avoids type mismatches (DiffEqNoiseProcess copies `TaskLocalRNG` into `Xoshiro`
   internally) and prevents duplicate random streams. For other RNG types,
   `seed = UInt64(0)` (sentinel for "no seed").
2. **`seed` provided (nonzero)**: Construct `Xoshiro(seed)`. `rng_provided = false`.
3. **Problem has a seed** (`prob.seed != 0`): Use it. `rng_provided = false`.
4. **Neither**: Generate a random seed and construct `Xoshiro` from it.
   `rng_provided = false`.

## Return values

- `rng::AbstractRNG` — the RNG to use for noise processes and `integrator.rng`.
- `seed::UInt64` — the seed for solution metadata (`RODESolution.seed::UInt64`).
  For `TaskLocalRNG` conversion, this is the derived seed. For other user-provided
  RNGs, `UInt64(0)` (no seed was used).
- `rng_provided::Bool` — `true` when the user explicitly passed an `rng` kwarg.
  Controls whether downstream consumers (JumpProcesses reseeding, noise process
  reseeding) should skip their own seed-based initialization.
"""
function _resolve_rng(rng, seed, prob)
    if rng !== nothing
        if !(rng isa Random.AbstractRNG)
            throw(
                ArgumentError(
                    "`rng` must be an `AbstractRNG`, got $(typeof(rng))."
                )
            )
        end
        # TaskLocalRNG is a zero-field singleton pointing to shared task-local
        # state. Storing it directly would cause type mismatches with noise
        # processes (DiffEqNoiseProcess converts it to Xoshiro internally) and
        # means set_rng! cannot sync integrator.rng with W.rng. Convert to a
        # concrete Xoshiro seeded from one draw of the task-local stream.
        if rng isa Random.TaskLocalRNG
            _seed = rand(rng, UInt64)
            return Random.Xoshiro(_seed), _seed, true
        end
        return rng, UInt64(0), true
    end
    _seed = if iszero(seed)
        if (!(prob isa DiffEqBase.AbstractRODEProblem) || iszero(prob.seed))
            seed_multiplier() * rand(UInt64)
        else
            prob.seed
        end
    else
        seed
    end
    return Random.Xoshiro(_seed), _seed, false
end

function DiffEqBase.solve!(integrator::SDEIntegrator)
    @inbounds while !isempty(integrator.opts.tstops)
        while integrator.tdir * integrator.t < first(integrator.opts.tstops)
            loopheader!(integrator)
            if integrator.do_error_check && check_error!(integrator) != ReturnCode.Success
                return integrator.sol
            end
            perform_step!(integrator, integrator.cache)
            loopfooter!(integrator)
            if isempty(integrator.opts.tstops)
                break
            end
        end
        handle_tstop!(integrator)
    end
    postamble!(integrator)

    f = integrator.sol.prob.f isa Tuple ? integrator.sol.prob.f[1] : integrator.sol.prob.f

    if DiffEqBase.has_analytic(f)
        DiffEqBase.calculate_solution_errors!(
            integrator.sol; timeseries_errors = integrator.opts.timeseries_errors,
            dense_errors = integrator.opts.dense_errors
        )
    end
    if integrator.sol.retcode != ReturnCode.Default
        return integrator.sol
    end
    return integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, ReturnCode.Success)
end

# Helpers

function handle_dt!(integrator)
    return if iszero(integrator.dt) && integrator.opts.adaptive
        auto_dt_reset!(integrator)
        if sign(integrator.dt) != integrator.tdir && !iszero(integrator.dt) &&
                !isnan(integrator.dt)
            error("Automatic dt setting has the wrong sign. Exiting. Please report this error.")
        end
        if isnan(integrator.dt)
            @SciMLMessage(
                "Automatic dt set the starting dt as NaN, causing instability.",
                integrator.opts.verbose, :dt_NaN
            )
        end
    elseif integrator.opts.adaptive && integrator.dt > zero(integrator.dt) &&
            integrator.tdir < 0
        integrator.dt *= integrator.tdir # Allow positive dt, but auto-convert
    end
end

function initialize_callbacks!(integrator, initialize_save = true)
    t = integrator.t
    u = integrator.u
    callbacks = integrator.opts.callback
    integrator.u_modified = true

    u_modified = initialize!(callbacks, u, t, integrator)

    # if the user modifies u, we need to fix previous values before initializing
    # FSAL in order for the starting derivatives to be correct
    if u_modified
        if isinplace(integrator.sol.prob)
            recursivecopy!(integrator.uprev, integrator.u)
        else
            integrator.uprev = integrator.u
        end

        if initialize_save &&
                (
                any((c) -> c.save_positions[2], callbacks.discrete_callbacks) ||
                    any((c) -> c.save_positions[2], callbacks.continuous_callbacks)
            )
            savevalues!(integrator, true)
        end
    end

    # reset this as it is now handled so the integrators should proceed as normal
    integrator.u_modified = false

    return if initialize_save
        SciMLBase.save_discretes_if_enabled!(integrator, integrator.opts.callback; skip_duplicates = true)
    end
end
