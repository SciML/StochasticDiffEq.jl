function solve(
  prob::AbstractRODEProblem,
alg::algType,timeseries=[],ts=[],ks=[],recompile::Type{Val{recompile_flag}}=Val{true};
kwargs...) where {algType<:Union{AbstractRODEAlgorithm,AbstractSDEAlgorithm},recompile_flag}

  integrator = init(prob,alg,timeseries,ts,ks,recompile;kwargs...)
  solve!(integrator)
  integrator.sol
end

function init(
  prob::AbstractRODEProblem{uType,tType,isinplace,ND},
  alg::algType,timeseries_init=uType[],ts_init=tType[],ks_init=[],
  recompile::Type{Val{recompile_flag}}=Val{true};
  dt = tType(0),
  timeseries_steps::Int = 1,
  save_noise = true,
  saveat = tType[],tstops = tType[],d_discontinuities= tType[],
  save_timeseries = nothing,
  save_everystep = isempty(saveat),
  save_idxs = nothing,
  save_start = true,save_end = true,
  dense = save_everystep,
  calck = (!isempty(setdiff(saveat,tstops)) || dense),
  adaptive=isadaptive(alg),gamma=9//10,
  abstol=1e-2,reltol=1e-2,
  qmax=qmax_default(alg),qmin=qmin_default(alg),
  failfactor = 2,
  qoldinit=1//10^4, fullnormalize=true,
  beta2=beta2_default(alg),
  beta1=beta1_default(alg,beta2),
  delta=1//6,maxiters = 1e9,
  dtmax=tType((prob.tspan[end]-prob.tspan[1])),
  dtmin=tType <: AbstractFloat ? tType(1000)*eps(tType) : tType(1//10^(10)),
  internalnorm=ODE_DEFAULT_NORM,
  unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
  isoutofdomain = ODE_DEFAULT_ISOUTOFDOMAIN,
  verbose = true,force_dtmin = false,
  advance_to_tstop = false,stop_at_next_tstop=false,
  progress_steps=1000,
  progress=false, progress_message = ODE_DEFAULT_PROG_MESSAGE,
  progress_name="SDE",
  userdata=nothing,callback=nothing,
  initialize_save = true,
  timeseries_errors = true, dense_errors=false,
  initialize_integrator=true,
  seed = UInt64(0),
  kwargs...) where {uType,tType,isinplace,algType<:Union{AbstractRODEAlgorithm,AbstractSDEAlgorithm},ND,recompile_flag}

  if save_timeseries != nothing
    warn("save_timeseries is deprecated. Use save_everystep instead")
    save_everystep = save_timeseries
  end

  if typeof(prob.f)<:Tuple
    if min((mm != I for mm in prob.mass_matrix)...)
      error("This solver is not able to use mass matrices.")
    end
  elseif prob.mass_matrix != I && !alg_mass_matrix_compatible(alg)
    error("This solver is not able to use mass matrices.")
  end

  if (typeof(prob.noise)<:NoiseProcess) && typeof(prob.noise.bridge)<:Void && adaptive
    error("Bridge function must be given for adaptivity. Either declare this function in noise process or set adaptive=false")
  end

  if !alg_compatible(prob,alg)
    error("The algorithm is not compatible with the chosen noise type. Please see the documentation on the solver methods")
  end

  noise = prob.noise
  tspan = prob.tspan
  tdir = sign(tspan[end]-tspan[1])

  T = tType(tspan[2])
  t = tType(tspan[1])

  if !adaptive && dt == tType(0) && isempty(tstops)
      error("Fixed timestep methods require a choice of dt or choosing the tstops")
  end

  f = prob.f
  if typeof(prob) <: AbstractSDEProblem
    g = prob.g
  else
    g = nothing
  end
  u0 = prob.u0
  uBottomEltype = recursive_bottom_eltype(u0)
  uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u0)

  if typeof(prob.u0) <: Array
    u = recursivecopy(prob.u0)
  elseif typeof(prob.u0) <: Number
    u = prob.u0
  elseif typeof(prob.u0) <: Tuple
    u = ArrayPartition(prob.u0,Val{true})
  else
    u = deepcopy(prob.u0)
  end

  ks = Vector{uType}(0)

  order = alg_order(alg)

  if typeof(alg) <: SRI || typeof(alg) <: SRA
    order = alg.tableau.order
  end

  dtmax > 0 && tdir < 0 && (dtmax *= tdir) # Allow positive dtmax, but auto-convert
  # dtmin is all abs => does not care about sign already.

  if typeof(u) <: AbstractArray
    rate_prototype = similar(u/zero(t),indices(u)) # rate doesn't need type info
  else
    rate_prototype = u/zero(t)
  end
  rateType = typeof(rate_prototype) ## Can be different if united

  if ND <: Void
    noise_rate_prototype = rate_prototype
  else
    if typeof(prob) <: AbstractSDEProblem
      noise_rate_prototype = prob.noise_rate_prototype
    else
      noise_rate_prototype = prob.rand_prototype
    end
  end

  tstops_internal, saveat_internal, d_discontinuities_internal =
    tstop_saveat_disc_handling(tstops,saveat,d_discontinuities,tdir,tspan,tType)

  callbacks_internal = CallbackSet(callback,prob.callback)

  uEltypeNoUnits = recursive_unitless_eltype(u0)
  tTypeNoUnits   = typeof(recursive_one(t))

  ### Algorithm-specific defaults ###
  # if save_idxs == nothing
  #   ksEltype = Vector{rateType}
  # else
  #   ks_prototype = rate_prototype[save_idxs]
  #   ksEltype = Vector{typeof(ks_prototype)}
  # end

  # Have to convert incase passed in wrong.
  if save_idxs == nothing
    timeseries = convert(Vector{uType},timeseries_init)
  else
    u_initial = u[save_idxs]
    timeseries = convert(Vector{typeof(u_initial)},timeseries_init)
  end
  ts = convert(Vector{tType},ts_init)
  #ks = convert(Vector{ksEltype},ks_init)
  alg_choice = Int[]

  if save_start
    saveiter = 1
    copyat_or_push!(ts,1,t)
    if save_idxs == nothing
      copyat_or_push!(timeseries,1,u)
      # copyat_or_push!(ks,1,[rate_prototype])
    else
      copyat_or_push!(timeseries,1,u_initial,Val{false})
      # copyat_or_push!(ks,1,[ks_prototype])
    end
  else
    saveiter = 0
  end

  opts = SDEOptions(maxiters,timeseries_steps,save_everystep,adaptive,map(uBottomEltype,abstol),
    map(uBottomEltypeNoUnits,reltol),tTypeNoUnits(gamma),tTypeNoUnits(qmax),tTypeNoUnits(qmin),
    tTypeNoUnits(failfactor),
    dtmax,dtmin,internalnorm,save_idxs,
    tstops_internal,saveat_internal,d_discontinuities_internal,
    tstops,saveat,d_discontinuities,
    userdata,
    progress,progress_steps,
    progress_name,progress_message,
    timeseries_errors,dense_errors,
    tTypeNoUnits(beta1),tTypeNoUnits(beta2),map(uBottomEltypeNoUnits,delta),tTypeNoUnits(qoldinit),
    dense,save_start,save_end,save_noise,
    callbacks_internal,isoutofdomain,unstable_check,verbose,calck,force_dtmin,
    advance_to_tstop,stop_at_next_tstop)

  progress ? (prog = Juno.ProgressBar(name=progress_name)) : prog = nothing

  # notsaveat_idxs = Int[1]

  # k = rateType[]

  if uType <: Array
    uprev = copy(u)
  else
    uprev = deepcopy(u)
  end

  dtcache = tType(dt)
  dtchangeable = true

  if !(uType <: AbstractArray)
    rand_prototype = zero(u/u) # Strip units and type info
    randType = typeof(rand_prototype)
  else
    randElType = typeof(recursive_one(u)) # Strip units and type info
    if ND <: Void # noise_dim isn't set, so it's diagonal
      if typeof(u) <: SArray
        rand_prototype = zero(u) # TODO: Array{randElType} for units
      else
        rand_prototype = similar(Array{randElType},indices(u))
        fill!(rand_prototype,zero(randElType))
      end
    elseif typeof(prob) <: AbstractSDEProblem
      rand_prototype = similar(Vector{randElType},size(noise_rate_prototype,2))
      fill!(rand_prototype,zero(randElType))
    else
      rand_prototype = prob.rand_prototype
    end
    randType = typeof(rand_prototype) # Strip units and type info
  end

  seed == 0 ? (prob.seed == 0 ? _seed = rand(UInt64) : _seed = prob.seed) : _seed = seed

  if typeof(prob.noise) <: Void
    if isinplace
      isadaptive(alg) ? rswm = RSWM(adaptivealg=:RSwM3) : rswm = RSWM(adaptivealg=:RSwM1)
      if alg_needs_extra_process(alg)
        W = WienerProcess!(t,rand_prototype,rand_prototype,
                           save_everystep=save_everystep,
                           timeseries_steps=timeseries_steps,
                           rswm=rswm,
                           rng = Xorshifts.Xoroshiro128Plus(_seed))
      else
        W = WienerProcess!(t,rand_prototype,
                           save_everystep=save_everystep,
                           timeseries_steps=timeseries_steps,
                           rswm=rswm,
                           rng = Xorshifts.Xoroshiro128Plus(_seed))
      end
    else
      if alg_needs_extra_process(alg)
        W = WienerProcess(t,rand_prototype,rand_prototype,
                           save_everystep=save_everystep,
                           timeseries_steps=timeseries_steps,
                           rng = Xorshifts.Xoroshiro128Plus(_seed))
      else
        W = WienerProcess(t,rand_prototype,
                           save_everystep=save_everystep,
                           timeseries_steps=timeseries_steps,
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
        srand(W.rng,_seed)
      end
    elseif W.t[end] != t
      error("Starting time in the noise process is not the starting time of the simulation. The noise process should be re-initialized for repeated use")
    end
  end

  EEst = tTypeNoUnits(1)
  q = tTypeNoUnits(1)
  just_hit_tstop = false
  isout = false
  accept_step = false
  u_modified = false
  last_stepfail = false
  force_stepfail = false
  tprev = t
  iter = 0
  success_iter = 0
  q11 = tTypeNoUnits(1)

  rateType = typeof(u/t) ## Can be different if united

  cache = alg_cache(alg,prob,u,W.dW,W.dZ,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,Val{isinplace})

  id = LinearInterpolationData(timeseries,ts)

  sol = build_solution(prob,alg,ts,timeseries,W=W,
                calculate_error = false,
                interp = id, dense = dense, seed = _seed)

  if recompile_flag == true
    FType = typeof(f)
    GType = typeof(g)
    SolType = typeof(sol)
    cacheType = typeof(cache)
  else
    FType = Function
    GType = Function
    SolType = AbstractODESolution
    cacheType =  OrdinaryDiffEqCache
  end

  integrator =    SDEIntegrator{typeof(alg),uType,uBottomEltype,tType,tTypeNoUnits,
                  uEltypeNoUnits,typeof(W),rateType,typeof(sol),typeof(cache),
                  typeof(prog),FType,GType,typeof(opts),typeof(noise)}(
                  f,g,noise,uprev,tprev,t,u,tType(dt),tType(dt),tType(dt),dtcache,T,tdir,
                  just_hit_tstop,isout,accept_step,last_stepfail,force_stepfail,
                  dtchangeable,u_modified,
                  saveiter,
                  alg,sol,
                  cache,tType(dt),W,
                  opts,iter,success_iter,prog,EEst,q,
                  tTypeNoUnits(qoldinit),q11)

  if initialize_integrator
    initialize_callbacks!(integrator)
    initialize!(integrator,integrator.cache)
  end

  if integrator.dt == zero(integrator.dt) && integrator.opts.adaptive
    auto_dt_reset!(integrator)
    if sign(integrator.dt)!=integrator.tdir && integrator.dt!=tType(0) && !isnan(integrator.dt)
      error("Automatic dt setting has the wrong sign. Exiting. Please report this error.")
    end
    if isnan(integrator.dt)
      if verbose
        warn("Automatic dt set the starting dt as NaN, causing instability.")
      end
    end
  elseif integrator.opts.adaptive && integrator.dt > zero(integrator.dt) && integrator.tdir < 0
    integrator.dt *= integrator.tdir # Allow positive dt, but auto-convert
  end

  ## Modify the first dt for tstops
  modify_dt_for_tstops!(integrator)
  ### Needs to be done before first rand
  integrator.sqdt = integrator.tdir*sqrt(abs(integrator.dt))

  integrator.W.dt = integrator.dt
  DiffEqNoiseProcess.setup_next_step!(integrator.W)


  integrator
end

function solve!(integrator::SDEIntegrator)
  @inbounds while !isempty(integrator.opts.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.opts.tstops)
      loopheader!(integrator)
      @sde_exit_condtions
      perform_step!(integrator,integrator.cache)
      loopfooter!(integrator)
      if isempty(integrator.opts.tstops)
        break
      end
    end
    handle_tstop!(integrator)
  end
  postamble!(integrator)

  if typeof(integrator.sol.prob.f) <: Tuple
    f = integrator.sol.prob.f[1]
  else
    f = integrator.sol.prob.f
  end

  if has_analytic(f)
    calculate_solution_errors!(integrator.sol;timeseries_errors=integrator.opts.timeseries_errors,dense_errors=integrator.opts.dense_errors)
  end
  integrator.sol = solution_new_retcode(integrator.sol,:Success)
  nothing
end

# Helpers

function tstop_saveat_disc_handling(tstops,saveat,d_discontinuities,tdir,tspan,tType)

  if isempty(d_discontinuities) && isempty(tstops) # TODO: Specialize more
    tstops_vec = [tspan[2]]
  else
    tstops_vec = vec(collect(tType,Iterators.filter(x->tdir*tspan[1]<tdir*x≤tdir*tspan[end],Iterators.flatten((tstops,d_discontinuities,tspan[end])))))
  end

  if tdir>0
    tstops_internal = binary_minheap(tstops_vec)
  else
    tstops_internal = binary_maxheap(tstops_vec)
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
    saveat_internal = binary_minheap(saveat_vec)
  else
    saveat_internal = binary_maxheap(saveat_vec)
  end

  d_discontinuities_vec = vec(collect(d_discontinuities))

  if tdir>0
    d_discontinuities_internal = binary_minheap(d_discontinuities_vec)
  else
    d_discontinuities_internal = binary_maxheap(d_discontinuities_vec)
  end
  tstops_internal,saveat_internal,d_discontinuities_internal
end

function initialize_callbacks!(integrator, initialize_save = true)
  t = integrator.t
  u = integrator.u
  callbacks = integrator.opts.callback
  integrator.u_modified = true

  u_modified = initialize!(callbacks,t,u,integrator)

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
