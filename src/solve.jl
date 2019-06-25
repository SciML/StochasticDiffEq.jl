function DiffEqBase.__solve(prob::DiffEqBase.AbstractRODEProblem,
                            alg::Union{AbstractRODEAlgorithm,AbstractSDEAlgorithm},
                            timeseries=[],ts=[],ks=nothing, # needed for variable rate
                            recompile::Type{Val{recompile_flag}}=Val{true};
                            kwargs...) where recompile_flag
  integrator = DiffEqBase.__init(prob,alg,timeseries,ts,recompile;kwargs...)
  solve!(integrator)
  integrator.sol
end

function DiffEqBase.__init(
  prob::DiffEqBase.AbstractRODEProblem,
  alg::Union{AbstractRODEAlgorithm,AbstractSDEAlgorithm},timeseries_init=typeof(prob.u0)[],
  ts_init=eltype(prob.tspan)[],
  ks_init=nothing,
  recompile::Type{Val{recompile_flag}}=Val{true};
  saveat = eltype(prob.tspan)[],
  tstops = eltype(prob.tspan)[],
  d_discontinuities= eltype(prob.tspan)[],
  save_idxs = nothing,
  save_everystep = isempty(saveat),
  save_noise = save_everystep && typeof(prob.f) <: Tuple ?
               DiffEqBase.has_analytic(prob.f[1]) : DiffEqBase.has_analytic(prob.f),
  save_on = true,
  save_start = save_everystep || isempty(saveat) || typeof(saveat) <: Number ? true : prob.tspan[1] in saveat,
  save_end = save_everystep || isempty(saveat) || typeof(saveat) <: Number ? true : prob.tspan[2] in saveat,
  callback=nothing,
  dense = save_everystep && isempty(saveat),
  calck = (!isempty(setdiff(saveat,tstops)) || dense),
  dt = eltype(prob.tspan)(0),
  adaptive = isadaptive(alg),
  gamma=9//10,
  abstol=nothing,
  reltol=nothing,
  qmax=qmax_default(alg),qmin=qmin_default(alg),
  qoldinit=1//10^4, fullnormalize=true,
  failfactor = 2,
  beta2=beta2_default(alg),
  beta1=beta1_default(alg,beta2),
  delta=delta_default(alg),
  maxiters = adaptive ? 1000000 : typemax(Int),
  dtmax=eltype(prob.tspan)((prob.tspan[end]-prob.tspan[1])),
  dtmin = typeof(one(eltype(prob.tspan))) <: AbstractFloat ? eps(eltype(prob.tspan)) :
          typeof(one(eltype(prob.tspan))) <: Integer ? 0 :
          eltype(prob.tspan)(1//10^(10)),
  internalnorm = ODE_DEFAULT_NORM,
  isoutofdomain = ODE_DEFAULT_ISOUTOFDOMAIN,
  unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
  verbose = true,force_dtmin = false,
  timeseries_errors = true, dense_errors=false,
  advance_to_tstop = false,stop_at_next_tstop=false,
  initialize_save = true,
  progress=false,progress_steps=1000,progress_name="SDE",
  progress_message = ODE_DEFAULT_PROG_MESSAGE,
  userdata=nothing,
  initialize_integrator=true,
  seed = UInt64(0), alias_u0=false, kwargs...) where recompile_flag

  if typeof(prob.f)<:Tuple
    if any(mm != I for mm in prob.f.mass_matrix)
      error("This solver is not able to use mass matrices.")
    end
  elseif prob.f.mass_matrix != I && !alg_mass_matrix_compatible(alg)
    error("This solver is not able to use mass matrices.")
  end

  if !isempty(saveat) && dense
    @warn("Dense output is incompatible with saveat. Please use the SavingCallback from the Callback Library to mix the two behaviors.")
  end

  if typeof(prob.noise)<:NoiseProcess && prob.noise.bridge === nothing && adaptive
    error("Bridge function must be given for adaptivity. Either declare this function in noise process or set adaptive=false")
  end

  if !alg_compatible(prob,alg)
    error("The algorithm is not compatible with the chosen noise type. Please see the documentation on the solver methods")
  end

  progress && @logmsg(-1,progress_name,_id=_id = :StochasticDiffEq,progress=0)

  tType = eltype(prob.tspan)
  noise = prob.noise
  tspan = prob.tspan
  tdir = sign(tspan[end]-tspan[1])

  t = tspan[1]

  if !adaptive && iszero(dt) && isempty(tstops)
      error("Fixed timestep methods require a choice of dt or choosing the tstops")
  end

  f = prob.f
  p = prob.p
  g = prob isa DiffEqBase.AbstractSDEProblem ? prob.g : nothing

  if typeof(prob.u0) <: Tuple
    u = ArrayPartition(prob.u0,Val{true})
  else
    if alias_u0
      u = prob.u0
    else
      u = recursivecopy(prob.u0)
    end
  end

  uType = typeof(u)
  uBottomEltype = recursive_bottom_eltype(u)
  uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)

  ks = Vector{uType}(undef, 0)

  uEltypeNoUnits = recursive_unitless_eltype(u)
  tTypeNoUnits   = typeof(one(tType))

  if abstol === nothing
    if uBottomEltypeNoUnits == uBottomEltype
      abstol_internal = real(convert(uBottomEltype,oneunit(uBottomEltype)*1//10^2))
    else
      abstol_internal = real.(oneunit.(u).*1//10^2)
    end
  else
    abstol_internal = real.(abstol)
  end

  if reltol === nothing
    if uBottomEltypeNoUnits == uBottomEltype
      reltol_internal = real(convert(uBottomEltype,oneunit(uBottomEltype)*1//10^2))
    else
      reltol_internal = real.(oneunit.(u).*1//10^2)
    end
  else
    reltol_internal = real.(reltol)
  end

  dtmax > zero(dtmax) && tdir < 0 && (dtmax *= tdir) # Allow positive dtmax, but auto-convert
  # dtmin is all abs => does not care about sign already.

  if isinplace(prob) && typeof(u) <: AbstractArray && eltype(u) <: Number && uBottomEltypeNoUnits == uBottomEltype # Could this be more efficient for other arrays?
    if !(typeof(u) <: ArrayPartition)
      rate_prototype = recursivecopy(u)
    else
      rate_prototype = similar(u, typeof.(oneunit.(recursive_bottom_eltype.(u.x))./oneunit(tType))...)
    end
  else
    if uBottomEltypeNoUnits == uBottomEltype
      rate_prototype = u
    else # has units!
      rate_prototype = u/oneunit(tType)
    end
  end
  rateType = typeof(rate_prototype) ## Can be different if united

  if is_diagonal_noise(prob)
    noise_rate_prototype = rate_prototype
  else
    if prob isa DiffEqBase.AbstractSDEProblem
      noise_rate_prototype = copy(prob.noise_rate_prototype)
    else
      noise_rate_prototype = copy(prob.rand_prototype)
    end
  end

  tstops_internal, saveat_internal, d_discontinuities_internal =
    tstop_saveat_disc_handling(tstops,saveat,d_discontinuities,tdir,tspan,tType)

  callbacks_internal = CallbackSet(callback,prob.callback)

  max_len_cb = DiffEqBase.max_vector_callback_length(callbacks_internal)
  if max_len_cb isa DiffEqBase.VectorContinuousCallback
    callback_cache = DiffEqBase.CallbackCache(max_len_cb.len,uBottomEltype,uBottomEltype)
  else
    callback_cache = nothing
  end

  ### Algorithm-specific defaults ###
  # if save_idxs === nothing
  #   ksEltype = Vector{rateType}
  # else
  #   ks_prototype = rate_prototype[save_idxs]
  #   ksEltype = Vector{typeof(ks_prototype)}
  # end

  # Have to convert incase passed in wrong.
  if save_idxs === nothing
    timeseries = convert(Vector{uType},timeseries_init)
  else
    u_initial = u[save_idxs]
    timeseries = convert(Vector{typeof(u_initial)},timeseries_init)
  end
  ts = convert(Vector{tType},ts_init)
  alg_choice = Int[]

  if !adaptive && save_everystep && tspan[2]-tspan[1] != Inf
    iszero(dt) ? steps = length(tstops) : steps = ceil(Int,internalnorm((tspan[2]-tspan[1])/dt,tspan[1]))
    sizehint!(timeseries,steps+1)
    sizehint!(ts,steps+1)
  elseif save_everystep
    sizehint!(timeseries,50)
    sizehint!(ts,50)
  elseif !isempty(saveat_internal)
    sizehint!(timeseries,length(saveat_internal)+1)
    sizehint!(ts,length(saveat_internal)+1)
  else
    sizehint!(timeseries,2)
    sizehint!(ts,2)
  end

  if save_start
    saveiter = 1 # Starts at 1 so first save is at 2
    copyat_or_push!(ts,1,t)
    if save_idxs === nothing
      copyat_or_push!(timeseries,1,u)
    else
      copyat_or_push!(timeseries,1,u_initial,Val{false})
    end
    if typeof(alg) <: StochasticDiffEqCompositeAlgorithm
      copyat_or_push!(alg_choice,1,1)
    end
  else
    saveiter = 0
  end

  QT = tTypeNoUnits <: Integer ? typeof(qmin) : tTypeNoUnits

  uprev = recursivecopy(u)

  if !(uType <: AbstractArray)
    rand_prototype = zero(u/u) # Strip units and type info
    randType = typeof(rand_prototype)
  else
    randElType = uBottomEltypeNoUnits # Strip units and type info
    if is_diagonal_noise(prob)
      if typeof(u) <: SArray
        rand_prototype = zero(u) # TODO: Array{randElType} for units
      else
        rand_prototype = (u .- u)./sqrt(oneunit(t))
      end
    elseif prob isa DiffEqBase.AbstractSDEProblem
      rand_prototype = noise_rate_prototype[1,:]
      fill!(rand_prototype,zero(randElType))
    else
      rand_prototype = copy(prob.rand_prototype)
    end
    randType = typeof(rand_prototype) # Strip units and type info
  end

  _seed = iszero(seed) ? (iszero(prob.seed) ? rand(UInt64) : prob.seed) : seed

  if prob.noise === nothing
    rswm = isadaptive(alg) ? RSWM(adaptivealg=:RSwM3) : RSWM(adaptivealg=:RSwM1)
    if isinplace(prob)
      if alg_needs_extra_process(alg)
        W = WienerProcess!(t,rand_prototype,rand_prototype,
                           save_everystep=save_noise,
                           rswm=rswm,
                           rng = Xorshifts.Xoroshiro128Plus(_seed))
      else
        W = WienerProcess!(t,rand_prototype,
                           save_everystep=save_noise,
                           rswm=rswm,
                           rng = Xorshifts.Xoroshiro128Plus(_seed))
      end
    else
      if alg_needs_extra_process(alg)
        W = WienerProcess(t,rand_prototype,rand_prototype,
                           save_everystep=save_noise,
                           rswm=rswm,
                           rng = Xorshifts.Xoroshiro128Plus(_seed))
      else
        W = WienerProcess(t,rand_prototype,
                           save_everystep=save_noise,
                           rswm=rswm,
                           rng = Xorshifts.Xoroshiro128Plus(_seed))
      end
    end
  else
    W = prob.noise
    if W.reset
      if W.t[end] != t
        reinit!(W,t)
      end
      # Reseed
      if typeof(W) <: NoiseProcess && W.reseed
        Random.seed!(W.rng,_seed)
      end
    elseif W.t[end] != t
      error("Starting time in the noise process is not the starting time of the simulation. The noise process should be re-initialized for repeated use")
    end
  end

  cache = alg_cache(alg,prob,u,W.dW,W.dZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,Val{isinplace(prob)})

  id = LinearInterpolationData(timeseries,ts)

  opts = SDEOptions(maxiters,save_everystep,
                    adaptive,abstol_internal,
                    reltol_internal,QT(gamma),
                    QT(qmax),QT(qmin),
                    QT(failfactor),
                    tType(dtmax),tType(dtmin),internalnorm,save_idxs,
                    tstops_internal,saveat_internal,
                    d_discontinuities_internal,
                    tstops,saveat,d_discontinuities,
                    userdata,
                    progress,progress_steps,
                    progress_name,progress_message,
                    timeseries_errors,dense_errors,
                    QT(beta1),QT(beta2),
                    convert.(uBottomEltypeNoUnits,delta),
                    QT(qoldinit),
                    dense,save_on,save_start,save_end,save_noise,
                    callbacks_internal,isoutofdomain,unstable_check,
                    verbose,calck,force_dtmin,
                    advance_to_tstop,stop_at_next_tstop)

  if typeof(alg) <: Union{StochasticDiffEqCompositeAlgorithm,
                          StochasticDiffEqRODECompositeAlgorithm}
    sol = DiffEqBase.build_solution(prob,alg,ts,timeseries,W=W,
                                    destats = DiffEqBase.DEStats(0),
                                    calculate_error = false, alg_choice=alg_choice,
                                    interp = id, dense = dense, seed = _seed)
  else
    sol = DiffEqBase.build_solution(prob,alg,ts,timeseries,W=W,
                                    destats = DiffEqBase.DEStats(0),
                                    calculate_error = false,
                                    interp = id, dense = dense, seed = _seed)
  end

  if recompile_flag == true
    FType = typeof(f)
    GType = typeof(g)
    SolType = typeof(sol)
    cacheType = typeof(cache)
  else
    FType = Function
    GType = Function
    SolType = DiffEqBase.AbstractRODESolution
    cacheType = StochasticDiffEqCache
  end

  tprev = t
  dtcache = tType(dt)
  iter = 0
  u_modified = false
  eigen_est = 1/oneunit(tType) # rate/state = (state/time)/state = 1/t units
  EEst = tTypeNoUnits(1)
  just_hit_tstop = false
  isout = false
  accept_step = false
  force_stepfail = false
  last_stepfail = false
  event_last_time = 0
  vector_event_last_time = 1
  last_event_error = zero(uBottomEltypeNoUnits)
  dtchangeable = true
  q11 = tTypeNoUnits(1)
  success_iter = 0
  q = tTypeNoUnits(1)

  integrator =    SDEIntegrator{typeof(alg),isinplace(prob),uType,
                  uBottomEltype,tType,typeof(p),
                  typeof(eigen_est),QT,
                  uEltypeNoUnits,typeof(W),rateType,typeof(sol),typeof(cache),
                  FType,GType,typeof(opts),typeof(noise),typeof(last_event_error),typeof(callback_cache)}(
                  f,g,noise,uprev,tprev,t,u,p,tType(dt),tType(dt),tType(dt),dtcache,tspan[2],tdir,
                  just_hit_tstop,isout,event_last_time,vector_event_last_time,last_event_error,accept_step,
                  last_stepfail,force_stepfail,
                  dtchangeable,u_modified,
                  saveiter,
                  alg,sol,
                  cache,callback_cache,tType(dt),W,
                  opts,iter,success_iter,eigen_est,EEst,q,
                  QT(qoldinit),q11)

  if initialize_integrator
    initialize_callbacks!(integrator, initialize_save)
    initialize!(integrator,integrator.cache)
    save_start && typeof(alg) <: Union{StochasticDiffEqCompositeAlgorithm,
                                       StochasticDiffEqRODECompositeAlgorithm} && copyat_or_push!(alg_choice,1,integrator.cache.current)
  end

  handle_dt!(integrator)

  ## Modify the first dt for tstops
  modify_dt_for_tstops!(integrator)
  ### Needs to be done before first rand
  integrator.sqdt = integrator.tdir*sqrt(abs(integrator.dt))

  integrator.W.dt = integrator.dt
  DiffEqNoiseProcess.setup_next_step!(integrator.W)

  integrator
end

function DiffEqBase.solve!(integrator::SDEIntegrator)
  @inbounds while !isempty(integrator.opts.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.opts.tstops)
      loopheader!(integrator)
      if check_error!(integrator) != :Success
        return integrator.sol
      end
      perform_step!(integrator,integrator.cache)
      loopfooter!(integrator)
      if isempty(integrator.opts.tstops)
        break
      end
    end
    handle_tstop!(integrator)
  end
  postamble!(integrator)

  f = typeof(integrator.sol.prob.f) <: Tuple ? integrator.sol.prob.f[1] : integrator.sol.prob.f

  if DiffEqBase.has_analytic(f)
    DiffEqBase.calculate_solution_errors!(integrator.sol;timeseries_errors=integrator.opts.timeseries_errors,dense_errors=integrator.opts.dense_errors)
  end
  if integrator.sol.retcode != :Default
    return integrator.sol
  end
  integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol,:Success)
end

# Helpers

function handle_dt!(integrator)
  if iszero(integrator.dt) && integrator.opts.adaptive
    auto_dt_reset!(integrator)
    if sign(integrator.dt)!=integrator.tdir && !iszero(integrator.dt) && !isnan(integrator.dt)
      error("Automatic dt setting has the wrong sign. Exiting. Please report this error.")
    end
    if isnan(integrator.dt)
      if integrator.opts.verbose
        @warn("Automatic dt set the starting dt as NaN, causing instability.")
      end
    end
  elseif integrator.opts.adaptive && integrator.dt > zero(integrator.dt) && integrator.tdir < 0
    integrator.dt *= integrator.tdir # Allow positive dt, but auto-convert
  end
end

function tstop_saveat_disc_handling(tstops,saveat,d_discontinuities,tdir,tspan,tType)

  if isempty(d_discontinuities) && isempty(tstops) # TODO: Specialize more
    tstops_vec = [tspan[2]]
  else
    tstops_vec = vec(collect(tType,Iterators.filter(x->tdir*tspan[1]<tdir*x≤tdir*tspan[end],Iterators.flatten((tstops,d_discontinuities,tspan[end])))))
  end

  if tdir>0
    tstops_internal = BinaryMinHeap(tstops_vec)
  else
    tstops_internal = BinaryMaxHeap(tstops_vec)
  end

  if typeof(saveat) <: Number
    if (tspan[1]:saveat:tspan[end])[end] == tspan[end]
      saveat_vec = convert(Vector{tType},collect(tType,tspan[1]+saveat:saveat:tspan[end]))
    else
      saveat_vec = convert(Vector{tType},collect(tType,tspan[1]+saveat:saveat:(tspan[end]-saveat)))
    end
  elseif isempty(saveat)
    saveat_vec = saveat
  else
    saveat_vec = vec(collect(tType,Iterators.filter(x->tdir*tspan[1]<tdir*x<tdir*tspan[end],saveat)))
  end

  if tdir>0
    saveat_internal = BinaryMinHeap(saveat_vec)
  else
    saveat_internal = BinaryMaxHeap(saveat_vec)
  end

  d_discontinuities_vec = vec(collect(d_discontinuities))

  if tdir>0
    d_discontinuities_internal = BinaryMinHeap(d_discontinuities_vec)
  else
    d_discontinuities_internal = BinaryMaxHeap(d_discontinuities_vec)
  end
  tstops_internal,saveat_internal,d_discontinuities_internal
end

function initialize_callbacks!(integrator, initialize_save = true)
  t = integrator.t
  u = integrator.u
  callbacks = integrator.opts.callback
  integrator.u_modified = true

  u_modified = initialize!(callbacks,u,t,integrator)

  # if the user modifies u, we need to fix previous values before initializing
  # FSAL in order for the starting derivatives to be correct
  if u_modified

    if isinplace(integrator.sol.prob)
      recursivecopy!(integrator.uprev,integrator.u)
    else
      integrator.uprev = integrator.u
    end

    if initialize_save &&
      (any((c)->c.save_positions[2],callbacks.discrete_callbacks) ||
      any((c)->c.save_positions[2],callbacks.continuous_callbacks))
      savevalues!(integrator,true)
    end
  end

  # reset this as it is now handled so the integrators should proceed as normal
  integrator.u_modified = false
end
