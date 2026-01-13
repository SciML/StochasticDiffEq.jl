@inline function DiffEqNoiseProcess.setup_next_step!(integrator::SDEIntegrator)
    !isnothing(integrator.W) &&
        DiffEqNoiseProcess.setup_next_step!(integrator.W, integrator.u, integrator.p)
    return !isnothing(integrator.P) &&
        DiffEqNoiseProcess.setup_next_step!(integrator.P, integrator.u, integrator.p)
end

@inline function DiffEqNoiseProcess.reject_step!(integrator::SDEIntegrator, dtnew = integrator.dtnew)
    !isnothing(integrator.W) &&
        reject_step!(integrator.W, dtnew, integrator.u, integrator.p)
    return !isnothing(integrator.P) &&
        reject_step!(integrator.P, dtnew, integrator.u, integrator.p)
end

@inline function DiffEqNoiseProcess.accept_step!(integrator::SDEIntegrator, setup)
    !isnothing(integrator.W) &&
        accept_step!(integrator.W, integrator.dt, integrator.u, integrator.p, setup)
    return !isnothing(integrator.P) &&
        accept_step!(integrator.P, integrator.dt, integrator.u, integrator.p, setup)
end

@inline function DiffEqNoiseProcess.save_noise!(integrator::SDEIntegrator)
    !isnothing(integrator.W) && DiffEqNoiseProcess.save_noise!(integrator.W)
    return !isnothing(integrator.P) && DiffEqNoiseProcess.save_noise!(integrator.P)
end

@inline function loopheader!(integrator::SDEIntegrator)
    # Apply right after iterators / callbacks

    # Accept or reject the step
    if integrator.iter > 0
        if (
                (integrator.opts.adaptive && integrator.accept_step) ||
                    !integrator.opts.adaptive || isaposteriori(integrator.alg)
            ) &&
                !integrator.force_stepfail
            integrator.success_iter += 1
            apply_step!(integrator)

        elseif integrator.opts.adaptive && !integrator.accept_step
            if integrator.isout
                integrator.dtnew = integrator.dt * integrator.opts.qmin
            elseif !integrator.force_stepfail
                step_reject_controller!(integrator, integrator.alg)
            end
            choose_algorithm!(integrator, integrator.cache)
            fix_dtnew_at_bounds!(integrator)
            modify_dtnew_for_tstops!(integrator)
            reject_step!(integrator)
            integrator.dt = integrator.dtnew
            integrator.sqdt = integrator.tdir * sqrt(abs(integrator.dt))
        end
    end

    integrator.iter += 1
    return integrator.force_stepfail = false
end

@inline function fix_dtnew_at_bounds!(integrator)
    integrator.dtnew = integrator.tdir * min(abs(integrator.opts.dtmax), abs(integrator.dtnew))
    return integrator.dtnew = integrator.tdir * max(abs(integrator.dtnew), abs(integrator.opts.dtmin))
end

@inline function modify_dt_for_tstops!(integrator)
    tstops = integrator.opts.tstops
    return @fastmath if !isempty(tstops)
        if integrator.opts.adaptive
            if integrator.tdir > 0
                integrator.dt = min(abs(integrator.dt), abs(first(tstops) - integrator.t)) # step! to the end
            else
                integrator.dt = -min(abs(integrator.dt), abs(first(tstops) + integrator.t))
            end
        elseif iszero(integrator.dtcache) && integrator.dtchangeable # Use integrator.opts.tstops
            integrator.dt = integrator.tdir *
                abs(first(tstops) - integrator.tdir * integrator.t)
        elseif integrator.dtchangeable && !integrator.force_stepfail
            # always try to step! with dtcache, but lower if a tstops
            integrator.dt = @fastmath integrator.tdir * min(
                abs(integrator.dtcache), abs(
                    first(tstops) -
                        integrator.tdir * integrator.t
                )
            ) # step! to the end
        end
    end
end

@inline function modify_dtnew_for_tstops!(integrator)
    tstops = integrator.opts.tstops
    return if !isempty(tstops)
        if integrator.tdir > 0
            integrator.dt = min(abs(integrator.dtnew), abs(first(tstops) - integrator.t)) # step! to the end
        else
            integrator.dt = -min(abs(integrator.dtnew), abs(first(tstops) + integrator.t))
        end
    end
end

function last_step_failed(integrator::SDEIntegrator)
    return integrator.last_stepfail && !integrator.opts.adaptive
end

@inline function savevalues!(integrator::SDEIntegrator, force_save = false)::Tuple{
        Bool, Bool,
    }
    saved, savedexactly = false, false
    !integrator.opts.save_on && return saved, savedexactly
    tdir_t = integrator.tdir * integrator.t
    while !isempty(integrator.opts.saveat) && first(integrator.opts.saveat) <= tdir_t # Perform saveat
        integrator.saveiter += 1
        saved = true
        curt = integrator.tdir * pop!(integrator.opts.saveat)
        if curt != integrator.t # If <t, interpolate
            Θ = (curt - integrator.tprev) / integrator.dt
            val = sde_interpolant(Θ, integrator, integrator.opts.save_idxs, Val{0}) # out of place, but force copy later
            save_val = val
            copyat_or_push!(integrator.sol.t, integrator.saveiter, curt)
            copyat_or_push!(integrator.sol.u, integrator.saveiter, save_val, Val{false})
            if integrator.alg isa StochasticDiffEqCompositeAlgorithm
                copyat_or_push!(integrator.sol.alg_choice, integrator.saveiter, integrator.cache.current)
            end
        else # ==t, just save
            savedexactly = true
            copyat_or_push!(integrator.sol.t, integrator.saveiter, integrator.t)
            if integrator.opts.save_idxs === nothing
                copyat_or_push!(integrator.sol.u, integrator.saveiter, integrator.u)
            else
                copyat_or_push!(
                    integrator.sol.u, integrator.saveiter,
                    integrator.u[integrator.opts.save_idxs], Val{false}
                )
            end
            if integrator.alg isa Union{
                    StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm,
                }
                copyat_or_push!(integrator.sol.alg_choice, integrator.saveiter, integrator.cache.current)
            end
        end
    end
    if force_save || (
            integrator.opts.save_everystep &&
                (isempty(integrator.sol.t) || (integrator.t !== integrator.sol.t[end]))
        )
        integrator.saveiter += 1
        saved, savedexactly = true, true
        if integrator.opts.save_idxs === nothing
            copyat_or_push!(integrator.sol.u, integrator.saveiter, integrator.u)
        else
            copyat_or_push!(
                integrator.sol.u, integrator.saveiter,
                integrator.u[integrator.opts.save_idxs], Val{false}
            )
        end
        copyat_or_push!(integrator.sol.t, integrator.saveiter, integrator.t)
        if integrator.alg isa
                Union{StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm}
            copyat_or_push!(integrator.sol.alg_choice, integrator.saveiter, integrator.cache.current)
        end
    end
    return saved, savedexactly
end

@inline function loopfooter!(integrator::SDEIntegrator)
    ttmp = integrator.t + integrator.dt
    integrator.do_error_check = true
    if integrator.force_stepfail
        if integrator.opts.adaptive
            integrator.dtnew = integrator.dt / integrator.opts.failfactor
        elseif integrator.last_stepfail
            return
        end
        integrator.last_stepfail = true
        integrator.accept_step = false
    elseif integrator.opts.adaptive
        stepsize_controller!(integrator, integrator.alg)
        integrator.isout = integrator.opts.isoutofdomain(integrator.u, integrator.p, ttmp)
        integrator.accept_step = (
            !integrator.isout &&
                accept_step_controller(integrator, integrator.opts.controller)
        ) ||
            (
            integrator.opts.force_dtmin &&
                integrator.dt <= integrator.opts.dtmin
        )
        if integrator.accept_step # Accepted
            step_accept_controller!(integrator, integrator.alg)
            integrator.last_stepfail = false
            integrator.tprev = integrator.t
            if typeof(integrator.t) <: AbstractFloat && !isempty(integrator.opts.tstops)
                tstop = integrator.tdir * first(integrator.opts.tstops)
                @fastmath abs(ttmp - tstop) < 100eps(integrator.t) ?
                    (integrator.t = tstop) : (integrator.t = ttmp)
            else
                integrator.t = ttmp
            end
            calc_dt_propose!(integrator)
            handle_callbacks!(integrator)
        end
    else # Non adaptive
        integrator.tprev = integrator.t
        if typeof(integrator.t) <: AbstractFloat && !isempty(integrator.opts.tstops)
            tstop = integrator.tdir * first(integrator.opts.tstops)
            # For some reason 100eps(integrator.t) is slow here
            # TODO: Allow higher precision but profile
            @fastmath abs(ttmp - tstop) < 100eps(max(integrator.t, tstop)) ?
                (integrator.t = tstop) : (integrator.t = ttmp)
        else
            integrator.t = ttmp
        end
        integrator.last_stepfail = false
        integrator.accept_step = true
        integrator.dtpropose = integrator.dt
        handle_callbacks!(integrator)
    end
    return if integrator.opts.progress && integrator.iter % integrator.opts.progress_steps == 0
        @logmsg(
            LogLevel(-1),
            integrator.opts.progress_name,
            _id = integrator.opts.progress_id,
            message = integrator.opts.progress_message(integrator.dt, integrator.u, integrator.p, integrator.t),
            progress = integrator.t / integrator.sol.prob.tspan[2]
        )
    end
end

@inline function calc_dt_propose!(integrator)
    integrator.qold = max(integrator.EEst, integrator.opts.qoldinit)
    if integrator.tdir > 0
        integrator.dtpropose = min(integrator.opts.dtmax, integrator.dtnew)
    else
        integrator.dtpropose = max(integrator.opts.dtmax, integrator.dtnew)
    end
    return if integrator.tdir > 0
        integrator.dtpropose = max(integrator.dtpropose, integrator.opts.dtmin) #abs to fix complex sqrt issue at end
    else
        integrator.dtpropose = min(integrator.dtpropose, integrator.opts.dtmin) #abs to fix complex sqrt issue at end
    end
end

@inline function solution_endpoint_match_cur_integrator!(integrator)
    if integrator.opts.save_end &&
            (
            integrator.saveiter == 0 ||
                integrator.sol.t[integrator.saveiter] != integrator.t &&
                (
                (integrator.opts.save_end_user isa Bool && integrator.opts.save_end_user) ||
                    integrator.t ∈ integrator.opts.saveat_cache ||
                    integrator.t == integrator.sol.prob.tspan[2] ||
                    isempty(integrator.opts.saveat_cache)
            )
        )
        integrator.saveiter += 1
        copyat_or_push!(integrator.sol.t, integrator.saveiter, integrator.t)
        if integrator.opts.save_idxs === nothing
            copyat_or_push!(integrator.sol.u, integrator.saveiter, integrator.u)
        else
            copyat_or_push!(
                integrator.sol.u, integrator.saveiter,
                integrator.u[integrator.opts.save_idxs], Val{false}
            )
        end
        SciMLBase.save_final_discretes!(integrator, integrator.opts.callback)
    end
    if (!isnothing(integrator.W) && integrator.W.curt != integrator.t) ||
            (!isnothing(integrator.P) && integrator.P.curt != integrator.t)
        accept_step!(integrator, false)
    end
    return if integrator.W isa NoiseProcess && !integrator.W.save_everystep
        save_noise!(integrator)
    end
end

@inline function postamble!(integrator::SDEIntegrator)
    DiffEqBase.finalize!(integrator.opts.callback, integrator.u, integrator.t, integrator)
    solution_endpoint_match_cur_integrator!(integrator)
    resize!(integrator.sol.t, integrator.saveiter)
    resize!(integrator.sol.u, integrator.saveiter)
    if integrator.opts.progress
        @logmsg(
            LogLevel(-1),
            integrator.opts.progress_name,
            _id = integrator.opts.progress_id,
            message = integrator.opts.progress_message(integrator.dt, integrator.u, integrator.p, integrator.t),
            progress = "done"
        )
    end
    return nothing
end

@inline function handle_callbacks!(integrator)
    discrete_callbacks = integrator.opts.callback.discrete_callbacks
    continuous_callbacks = integrator.opts.callback.continuous_callbacks
    atleast_one_callback = false

    continuous_modified = false
    discrete_modified = false
    saved_in_cb = false
    if !(typeof(continuous_callbacks) <: Tuple{})
        time, upcrossing,
            event_occurred,
            event_idx,
            idx,
            counter = DiffEqBase.find_first_continuous_callback(integrator, continuous_callbacks...)
        if event_occurred
            integrator.event_last_time = idx
            integrator.vector_event_last_time = event_idx
            continuous_modified,
                saved_in_cb = DiffEqBase.apply_callback!(
                integrator, continuous_callbacks[idx], time, upcrossing, event_idx
            )
        else
            integrator.event_last_time = 0
            integrator.vector_event_last_time = 1
        end
    end
    if !(typeof(discrete_callbacks) <: Tuple{})
        discrete_modified,
            saved_in_cb = DiffEqBase.apply_discrete_callback!(integrator, discrete_callbacks...)
    end
    if !saved_in_cb
        savevalues!(integrator)
    end

    integrator.u_modified = continuous_modified || discrete_modified
    return if integrator.u_modified
        integrator.do_error_check = false
        handle_callback_modifiers!(integrator)
    end
end

@inline function handle_callback_modifiers!(integrator::SDEIntegrator)
    #integrator.reeval_fsal = true
    return if integrator.P !== nothing && integrator.opts.adaptive
        if integrator.cache isa StochasticDiffEqMutableCache
            oldrate = integrator.P.cache.currate
            integrator.P.cache.rate(oldrate, integrator.u, integrator.p, integrator.t)
        else
            integrator.P.cache.currate = integrator.P.cache.rate(integrator.u, integrator.p, integrator.t)
        end
    end
end

@inline function apply_step!(integrator)
    if isinplace(integrator.sol.prob)
        recursivecopy!(integrator.uprev, integrator.u)
    else
        integrator.uprev = integrator.u
    end
    integrator.dt = integrator.dtpropose
    modify_dt_for_tstops!(integrator)
    accept_step!(integrator, true)

    # Allow RSWM1 on Wiener Process to change dt
    !isnothing(integrator.W) && (integrator.dt = integrator.W.dt)
    return integrator.sqdt = @fastmath integrator.tdir * sqrt(abs(integrator.dt)) # It can change dt, like in RSwM1
end

@inline function handle_tstop!(integrator)
    tstops = integrator.opts.tstops
    return if !isempty(tstops)
        tdir_t = integrator.tdir * integrator.t
        tdir_ts_top = first(tstops)
        if tdir_t == tdir_ts_top
            pop!(tstops)
            integrator.just_hit_tstop = true
        elseif tdir_t > tdir_ts_top
            if !integrator.dtchangeable
                change_t_via_interpolation!(integrator, integrator.tdir * pop!(tstops), Val{true})
                integrator.just_hit_tstop = true
            else
                error("Something went wrong. Integrator stepped past tstops but the algorithm was dtchangeable. Please report this error.")
            end
        end
    end
end

@inline function update_noise!(integrator, scaling_factor = integrator.sqdt)
    return if isinplace(integrator.noise)
        integrator.noise(integrator.ΔW, integrator)
        rmul!(integrator.ΔW, scaling_factor)
        if alg_needs_extra_process(integrator.alg)
            integrator.noise(integrator.ΔZ, integrator)
            rmul!(integrator.ΔZ, scaling_factor)
        end
    else
        if integrator.u isa AbstractArray
            integrator.ΔW .= scaling_factor .* integrator.noise(size(integrator.u), integrator)
            if alg_needs_extra_process(integrator.alg)
                integrator.ΔZ .= scaling_factor .* integrator.noise(size(integrator.u), integrator)
            end
        else
            integrator.ΔW = scaling_factor * integrator.noise(integrator)
            if alg_needs_extra_process(integrator.alg)
                integrator.ΔZ = scaling_factor * integrator.noise(integrator)
            end
        end
    end
end

@inline function generate_tildes(integrator, add1, add2, scaling)
    return if isinplace(integrator.noise)
        integrator.noise(integrator.ΔWtilde, integrator)
        if add1 != 0
            @.. integrator.ΔWtilde = add1 + scaling * integrator.ΔWtilde
        else
            @.. integrator.ΔWtilde = scaling * integrator.ΔWtilde
        end
        if alg_needs_extra_process(integrator.alg)
            integrator.noise(integrator.ΔZtilde, integrator)
            if add2 != 0
                @.. integrator.ΔZtilde = add2 + scaling * integrator.ΔZtilde
            else
                @.. integrator.ΔZtilde = scaling * integrator.ΔZtilde
            end
        end
    else
        if integrator.u isa AbstractArray
            if add1 != 0
                integrator.ΔWtilde = add1 .+ scaling .* integrator.noise(size(integrator.u), integrator)
            else
                integrator.ΔWtilde = scaling .* integrator.noise(size(integrator.u), integrator)
            end
            if alg_needs_extra_process(integrator.alg)
                if add2 != 0
                    integrator.ΔZtilde = add2 .+ scaling .* integrator.noise(size(integrator.u), integrator)
                else
                    integrator.ΔZtilde = scaling .* integrator.noise(size(integrator.u), integrator)
                end
            end
        else
            integrator.ΔWtilde = add1 + scaling * integrator.noise(integrator)
            if alg_needs_extra_process(integrator.alg)
                integrator.ΔZtilde = add2 + scaling * integrator.noise(integrator)
            end
        end
    end
end

@inline initialize!(integrator, cache::StochasticDiffEqCache, f = integrator.f) = nothing

function nlsolve!(integrator, cache)
    return DiffEqBase.nlsolve!(cache.nlsolver, cache.nlsolver.cache, integrator)
end

function OrdinaryDiffEqCore.nlsolve_f(f, alg::StochasticDiffEqAlgorithm)
    return f isa SplitSDEFunction && issplit(alg) ? f.f1 : f
end
function OrdinaryDiffEqCore.nlsolve_f(integrator::SDEIntegrator)
    return nlsolve_f(integrator.f, unwrap_alg(integrator, true))
end

# TauLeapingDrift: wrapper for tau-leaping drift function used by nlsolver
# Computes drift(u, p, t) = c(u, p, t, rate(u, p, t), nothing)
# where c is the stoichiometry function and rate is the propensity function
struct TauLeapingDrift{C, R, RateCache, IIP}
    c::C              # Stoichiometry function (from integrator.c)
    rate::R           # Rate function (from integrator.P.cache.rate)
    rate_cache::RateCache  # Cache for rate values (for in-place version)
end

# Out-of-place: drift(u, p, t)
function (td::TauLeapingDrift{C, R, Nothing, false})(u, p, t) where {C, R}
    rates = td.rate(u, p, t)
    return td.c(u, p, t, rates, nothing)
end

# In-place: drift(du, u, p, t)
function (td::TauLeapingDrift{C, R, RateCache, true})(du, u, p, t) where {C, R, RateCache}
    td.rate(td.rate_cache, u, p, t)
    td.c(du, u, p, t, td.rate_cache, nothing)
    return nothing
end

# nlsolve_f override for ImplicitTauLeaping: return TauLeapingDrift wrapper
function OrdinaryDiffEqCore.nlsolve_f(integrator::SDEIntegrator{A}) where {A <: ImplicitTauLeaping}
    # Determine if the cache is in-place or constant (out-of-place) based on cache type
    cache = integrator.cache
    if cache isa ImplicitTauLeapingCache
        # In-place version - use rate cache from the integrator's cache
        rate_cache = cache.rate_at_uprev
        return TauLeapingDrift{typeof(integrator.c), typeof(integrator.P.cache.rate), typeof(rate_cache), true}(
            integrator.c, integrator.P.cache.rate, rate_cache
        )
    else
        # Out-of-place (constant cache) version
        return TauLeapingDrift{typeof(integrator.c), typeof(integrator.P.cache.rate), Nothing, false}(
            integrator.c, integrator.P.cache.rate, nothing
        )
    end
end

# nlsolve_f override for ThetaTrapezoidalTauLeaping: return TauLeapingDrift wrapper
function OrdinaryDiffEqCore.nlsolve_f(integrator::SDEIntegrator{A}) where {A <: ThetaTrapezoidalTauLeaping}
    # Determine if the cache is in-place or constant (out-of-place) based on cache type
    cache = integrator.cache
    if cache isa ThetaTrapezoidalTauLeapingCache
        # In-place version - use rate cache from the integrator's cache
        rate_cache = cache.rate_at_uprev
        return TauLeapingDrift{typeof(integrator.c), typeof(integrator.P.cache.rate), typeof(rate_cache), true}(
            integrator.c, integrator.P.cache.rate, rate_cache
        )
    else
        # Out-of-place (constant cache) version
        return TauLeapingDrift{typeof(integrator.c), typeof(integrator.P.cache.rate), Nothing, false}(
            integrator.c, integrator.P.cache.rate, nothing
        )
    end
end

function iip_generate_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits)
    if alg.nlsolve isa NLNewton
        nf = nlsolve_f(f, alg)
        islin = f isa Union{SDEFunction, SplitSDEFunction} && islinear(nf.f)
        if islin
            J = nf.f
            W = WOperator{true}(f.mass_matrix, dt, J, u)
        else
            if ArrayInterface.isstructured(f.jac_prototype) ||
                    f.jac_prototype isa SparseMatrixCSC
                J = similar(f.jac_prototype)
                W = similar(J)
            elseif DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) &&
                    f.jac_prototype !== nothing
                J = nothing
                W = WOperator{true}(f, u, dt)
            else
                J = false .* vec(u) .* vec(u)'
                W = similar(J)
            end
        end
    else
        J = nothing
        W = nothing
    end
    return J, W
end

function oop_generate_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits)
    nf = nlsolve_f(f, alg)
    islin = f isa Union{SDEFunction, SplitSDEFunction} && islinear(nf.f)
    if islin || DiffEqBase.has_jac(f)
        # get the operator
        J = islin ? nf.f : f.jac(uprev, p, t)
        if !isa(J, DiffEqBase.AbstractDiffEqLinearOperator)
            J = MatrixOperator(J)
        end
        W = WOperator{false}(f.mass_matrix, dt, J, u)
    else
        if u isa StaticArray
            # get a "fake" `J`
            J = if u isa AbstractMatrix && size(u, 1) > 1 # `u` is already a matrix
                u
            elseif size(u, 1) == 1 # `u` is a row vector
                vcat(u, u)
            else # `u` is a column vector
                hcat(u, u)
            end
            W = lu(J)
        else
            W = u isa Number ? u :
                LU{LinearAlgebra.lutype(uEltypeNoUnits)}(
                    Matrix{uEltypeNoUnits}(undef, 0, 0),
                    Vector{LinearAlgebra.BlasInt}(undef, 0),
                    zero(LinearAlgebra.BlasInt)
                )
            J = u isa Number ? u : (false .* vec(u) .* vec(u)')
        end
    end
    return J, W
end
