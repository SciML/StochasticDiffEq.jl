@inline function DiffEqBase.change_t_via_interpolation!(integrator::SDEIntegrator,t,modify_save_endpoint::Type{Val{T}}=Val{false}) where T
  # Can get rid of an allocation here with a function
  # get_tmp_arr(integrator.cache) which gives a pointer to some
  # cache array which can be modified.
  if integrator.tdir*t < integrator.tdir*integrator.tprev
    error("Current interpolant only works between tprev and t")
  elseif t != integrator.t
    if typeof(integrator.u) <: AbstractArray
      integrator(integrator.u,t)
    else
      integrator.u = integrator(t)
    end
    integrator.dtnew = integrator.t - t
    reject_step!(integrator.W,t-integrator.tprev) #this only changes dt and noise, so no interpolation problems
    integrator.dt = integrator.dtnew
    integrator.sqdt = sqrt(abs(integrator.dt))
    integrator.t = t
    # reeval_internals_due_to_modification!(integrator) # Not necessary for linear interp
    if T
      solution_endpoint_match_cur_integrator!(integrator)
    end
  end
end

function (integrator::SDEIntegrator)(t,deriv::Type=Val{0};idxs=nothing)
  current_interpolant(t,integrator,idxs,deriv)
end

(integrator::SDEIntegrator)(val::AbstractArray,t::Union{Number,AbstractArray},deriv::Type=Val{0};idxs=nothing) = current_interpolant!(val,t,integrator,idxs,deriv)

function u_modified!(integrator::SDEIntegrator,bool::Bool)
  integrator.u_modified = bool
end

get_proposed_dt(integrator::SDEIntegrator) = integrator.dtpropose
set_proposed_dt!(integrator::SDEIntegrator,dt::Number) = (integrator.dtpropose = dt)

function set_proposed_dt!(integrator::SDEIntegrator,integrator2::SDEIntegrator)
  integrator.dtpropose = integrator2.dtpropose
  integrator.qold = integrator2.qold
  integrator.erracc = integrator2.erracc
  integrator.dtacc = integrator2.dtacc
end

#TODO: Bigger caches for most algorithms
@inline DiffEqBase.get_tmp_cache(integrator::SDEIntegrator) =
  get_tmp_cache(integrator, integrator.alg, integrator.cache)
# avoid method ambiguity
for typ in (StochasticDiffEqAlgorithm,StochasticDiffEqNewtonAdaptiveAlgorithm)
  @eval @inline DiffEqBase.get_tmp_cache(integrator::SDEIntegrator, alg::$typ, cache::StochasticDiffEqConstantCache) = nothing
end
@inline DiffEqBase.get_tmp_cache(integrator::SDEIntegrator, alg, cache) = (cache.tmp,)
@inline DiffEqBase.get_tmp_cache(integrator::SDEIntegrator, alg::StochasticDiffEqNewtonAdaptiveAlgorithm, cache) =
    (cache.tmp, cache.atmp)
@inline DiffEqBase.get_tmp_cache(integrator::SDEIntegrator, alg::StochasticCompositeAlgorithm, cache) =
    get_tmp_cache(integrator, alg.algs[1], cache.caches[1])

full_cache(integrator::SDEIntegrator) = full_cache(integrator.cache)
ratenoise_cache(integrator::SDEIntegrator) = ratenoise_cache(integrator.cache)
rand_cache(integrator::SDEIntegrator) = rand_cache(integrator.cache)
jac_iter(integrator::SDEIntegrator) = jac_iter(integrator.cache)

@inline function add_tstop!(integrator::SDEIntegrator,t)
  t < integrator.t && error("Tried to add a tstop that is behind the current time. This is strictly forbidden")
  push!(integrator.opts.tstops,t)
end

function DiffEqBase.add_saveat!(integrator::SDEIntegrator,t)
  integrator.tdir * (t - integrator.t) < 0 && error("Tried to add a saveat that is behind the current time. This is strictly forbidden")
  push!(integrator.opts.saveat,t)
end

resize_non_user_cache!(integrator::SDEIntegrator,i::Int) = resize_non_user_cache!(integrator,integrator.cache,i)
deleteat_non_user_cache!(integrator::SDEIntegrator,i) = deleteat_non_user_cache!(integrator,integrator.cache,i)
addat_non_user_cache!(integrator::SDEIntegrator,i) = addat_non_user_cache!(integrator,integrator.cache,i)
resize!(integrator::SDEIntegrator,i::Int) = resize!(integrator,integrator.cache,i)

function resize!(integrator::SDEIntegrator,cache,i)
  # This has to go first!
  resize_non_user_cache!(integrator,cache,i)
  for c in full_cache(integrator)
    resize!(c,i)
  end
  for c in ratenoise_cache(integrator)
    resize!(c,i)
  end
end

function resize_noise!(integrator,cache,bot_idx,i)
  for c in integrator.W.S₁
    resize!(c[2],i)
    if alg_needs_extra_process(integrator.alg)
      resize!(c[3],i)
    end
    if i >= bot_idx # fill in rands
      fill_new_noise_caches!(integrator,c,c[1],bot_idx:i)
    end
  end
  for c in integrator.W.S₂
    resize!(c[2],i)
    if alg_needs_extra_process(integrator.alg)
      resize!(c[3],i)
    end
    if i >= bot_idx # fill in rands
      fill_new_noise_caches!(integrator,c,c[1],bot_idx:i)
    end
  end
  resize!(integrator.W.dW,i)
  integrator.W.dW[end] = zero(eltype(integrator.u))
  resize!(integrator.W.dWtilde,i)
  integrator.W.dWtilde[end] = zero(eltype(integrator.u))
  resize!(integrator.W.dWtmp,i)
  integrator.W.dWtmp[end] = zero(eltype(integrator.u))
  resize!(integrator.W.curW,i)
  integrator.W.curW[end] = zero(eltype(integrator.u))
  DiffEqNoiseProcess.resize_stack!(integrator.W,i)

  if alg_needs_extra_process(integrator.alg)
    resize!(integrator.W.dZ,i)
    integrator.W.dZ[end] = zero(eltype(integrator.u))
    resize!(integrator.W.dZtilde,i)
    integrator.W.dZtilde[end] = zero(eltype(integrator.u))
    resize!(integrator.W.dZtmp,i)
    integrator.W.dZtmp[end] = zero(eltype(integrator.u))
    resize!(integrator.W.curZ,i)
    integrator.W.curZ[end] = zero(eltype(integrator.u))
  end
  if i >= bot_idx # fill in rands
    fill!(@view(integrator.W.curW[bot_idx:i]),zero(eltype(integrator.u)))
    if alg_needs_extra_process(integrator.alg)
      fill!(@view(integrator.W.curZ[bot_idx:i]),zero(eltype(integrator.u)))
    end
  end
end

@inline function fill_new_noise_caches!(integrator,c,scaling_factor,idxs)
  if isinplace(integrator.W)
    integrator.W.dist(@view(c[2][idxs]),integrator.W,scaling_factor,integrator.W.rng)
    if alg_needs_extra_process(integrator.alg)
      integrator.W.dist(@view(c[3][idxs]),integrator.W,scaling_factor,integrator.W.rng)
    end
  else
    c[2][idxs] .= integrator.noise(length(idxs),integrator,scaling_factor)
    if alg_needs_extra_process(integrator.alg)
      c[3][idxs] .= integrator.noise(length(idxs),integrator,scaling_factor)
    end
  end

end

function resize_non_user_cache!(integrator::SDEIntegrator,cache,i)
  bot_idx = length(integrator.u) + 1
  if is_diagonal_noise(integrator.sol.prob)
    resize_noise!(integrator,cache,bot_idx,i)
    for c in rand_cache(integrator)
      resize!(c,i)
    end
  end
end

function deleteat!(integrator::SDEIntegrator,idxs)
  deleteat_non_user_cache!(integrator,cache,idxs)
  for c in full_cache(integrator)
    deleteat!(c,idxs)
  end
  for c in ratenoise_cache(integrator)
    deleteat!(c,idxs)
  end
end

function addat!(integrator::SDEIntegrator,idxs)
  addat_non_user_cache!(integrator,cache,idxs)
  for c in full_cache(integrator)
    addat!(c,idxs)
  end
  for c in ratenoise_cache(integrator)
    addat!(c,idxs)
  end
end

function deleteat_non_user_cache!(integrator::SDEIntegrator,cache,idxs)
  if is_diagonal_noise(integrator.sol.prob)
    deleteat_noise!(integrator,cache,idxs)
    for c in rand_cache(integrator)
      deleteat!(c,idxs)
    end
  end
end

function addat_non_user_cache!(integrator::SDEIntegrator,cache,idxs)
  if is_diagonal_noise(integrator.sol.prob)
    addat_noise!(integrator,cache,idxs)
    for c in rand_cache(integrator)
      addat!(c,idxs)
    end
  end
end

function deleteat_noise!(integrator,cache,idxs)
  for c in integrator.W.S₁
    deleteat!(c[2],idxs)
    if alg_needs_extra_process(integrator.alg)
      deleteat!(c[3],idxs)
    end
  end
  for c in integrator.W.S₂
    deleteat!(c[2],idxs)
    if alg_needs_extra_process(integrator.alg)
      deleteat!(c[3],idxs)
    end
  end
  deleteat!(integrator.W.dW,idxs)
  deleteat!(integrator.W.dWtilde,idxs)
  deleteat!(integrator.W.dWtmp,idxs)
  deleteat!(integrator.W.curW,idxs)
  DiffEqNoiseProcess.resize_stack!(integrator.W,length(integrator.u))

  if alg_needs_extra_process(integrator.alg)
    deleteat!(integrator.W.curZ,idxs)
    deleteat!(integrator.W.dZtmp,idxs)
    deleteat!(integrator.W.dZtilde,idxs)
    deleteat!(integrator.W.dZ,idxs)
  end
end

function addat_noise!(integrator,cache,idxs)
  for c in integrator.W.S₁
    addat!(c[2],idxs)
    if alg_needs_extra_process(integrator.alg)
      addat!(c[3],idxs)
    end
    fill_new_noise_caches!(integrator,c,c[1],idxs)
  end
  for c in integrator.W.S₂
    addat!(c[2],idxs)
    if alg_needs_extra_process(integrator.alg)
      addat!(c[3],idxs)
    end
    fill_new_noise_caches!(integrator,c,c[1],idxs)
  end

  addat!(integrator.W.dW,idxs)
  integrator.W.dW[idxs] .= zero(eltype(integrator.u))
  addat!(integrator.W.curW,idxs)
  integrator.W.curW[idxs] .= zero(eltype(integrator.u))
  if alg_needs_extra_process(integrator.alg)
    addat!(integrator.W.dZ,idxs)
    integrator.W.dZ[idxs] .= zero(eltype(integrator.u))
    addat!(integrator.W.curZ,idxs)
    integrator.W.curZ[idxs] .= zero(eltype(integrator.u))
  end

  i = length(integrator.u)
  resize!(integrator.W.dWtilde,i)
  resize!(integrator.W.dWtmp,i)
  DiffEqNoiseProcess.resize_stack!(integrator.W,i)
  if alg_needs_extra_process(integrator.alg)
    resize!(integrator.W.dZtmp,i)
    resize!(integrator.W.dZtilde,i)
  end

  # fill in rands
  fill!(@view(integrator.W.curW[idxs]),zero(eltype(integrator.u)))
  if alg_needs_extra_process(integrator.alg)
    fill!(@view(integrator.W.curZ[idxs]),zero(eltype(integrator.u)))
  end
end


function terminate!(integrator::SDEIntegrator, retcode = :Terminated)
  integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, retcode)
  integrator.opts.tstops.valtree = typeof(integrator.opts.tstops.valtree)()
end

DiffEqBase.has_reinit(integrator::SDEIntegrator) = true
function DiffEqBase.reinit!(integrator::SDEIntegrator,u0 = integrator.sol.prob.u0;
  t0 = integrator.sol.prob.tspan[1], tf = integrator.sol.prob.tspan[2],
  erase_sol = true,
  tstops = integrator.opts.tstops_cache,
  saveat = integrator.opts.saveat_cache,
  d_discontinuities = integrator.opts.d_discontinuities_cache,
  reinit_cache = true,reinit_callbacks = true,
  initialize_save = true,
  reset_dt = (integrator.dtcache == zero(integrator.dt)) && integrator.opts.adaptive)

  if isinplace(integrator.sol.prob)
    recursivecopy!(integrator.u,u0)
    recursivecopy!(integrator.uprev,integrator.u)
  else
    integrator.u = u0
    integrator.uprev = integrator.u
  end

  integrator.t = t0
  integrator.tprev = t0

  tstops_internal, saveat_internal, d_discontinuities_internal =
    tstop_saveat_disc_handling(tstops,saveat,d_discontinuities,
    integrator.tdir,(t0,tf),typeof(integrator.t))

  integrator.opts.tstops = tstops_internal
  integrator.opts.saveat = saveat_internal
  integrator.opts.d_discontinuities = d_discontinuities_internal

  if erase_sol
    if integrator.opts.save_start
      resize_start = 1
    else
      resize_start = 0
    end
    resize!(integrator.sol.u,resize_start)
    resize!(integrator.sol.t,resize_start)
    if integrator.sol.u_analytic !== nothing
      resize!(integrator.sol.u_analytic,0)
    end
    if typeof(integrator.alg) <: StochasticDiffEqCompositeAlgorithm
      resize!(integrator.sol.alg_choice,resize_start)
    end
    integrator.saveiter = resize_start
  end
  integrator.iter = 0
  integrator.success_iter = 0

  # full re-initialize the PI in timestepping
  integrator.qold = integrator.opts.qoldinit
  integrator.q11 = typeof(integrator.t)(1)

  if reset_dt
    auto_dt_reset!(integrator)
  end

  if reinit_callbacks
    initialize_callbacks!(integrator, initialize_save)
  end

  if reinit_cache
    initialize!(integrator,integrator.cache)
  end

  reinit!(integrator.W,integrator.dt)
end

function DiffEqBase.auto_dt_reset!(integrator::SDEIntegrator)
  integrator.dt = sde_determine_initdt(integrator.u,integrator.t,
  integrator.tdir,integrator.opts.dtmax,integrator.opts.abstol,integrator.opts.reltol,
  integrator.opts.internalnorm,integrator.sol.prob,get_current_alg_order(integrator.alg, integrator.cache),
  integrator)
end

@inline function DiffEqBase.get_du(integrator::SDEIntegrator)
  (integrator.u - integrator.uprev) / integrator.dt
end

@inline function DiffEqBase.get_du!(out,integrator::SDEIntegrator)
  @.. out = (integrator.u - integrator.uprev) / integrator.dt
end

function DiffEqBase.set_t!(integrator::SDEIntegrator, t::Real)
  if integrator.opts.save_everystep
    error("Integrator time cannot be reset unless it is initialized",
          " with save_everystep=false")
  end
  if !isdtchangeable(integrator.alg)
    reinit!(integrator, integrator.u;
            t0 = t,
            reset_dt = false,
            reinit_callbacks = false,
            reinit_cache = false)
  else
    integrator.t = t
  end
end

function DiffEqBase.set_u!(integrator::SDEIntegrator, u)
  if integrator.opts.save_everystep
    error("Integrator state cannot be reset unless it is initialized",
          " with save_everystep=false")
  end
  integrator.u = u
  u_modified!(integrator, true)
end
