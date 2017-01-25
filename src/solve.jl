@inline ODE_DEFAULT_NORM(u) = sqrt(sum(abs2,u) / length(u))
@inline ODE_DEFAULT_PROG_MESSAGE(dt,t,u) = "dt="*string(dt)*"\nt="*string(t)*"\nmax u="*string(maximum(abs.(u)))
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,t,u) = any(isnan,u)

function solve{uType,tType,isinplace,NoiseClass,F,F2,F3,algType<:AbstractSDEAlgorithm,recompile_flag}(
              prob::AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3},
              alg::algType,timeseries_init=uType[],ts_init=tType[],ks_init=[],
              recompile::Type{Val{recompile_flag}}=Val{true};
              dt = tType(0),save_timeseries::Bool = true,
              timeseries_steps::Int = 1,
              dense = false,
              save_noise = true,
              saveat = tType[],tstops = tType[],d_discontinuities= tType[],
              calck = (!isempty(setdiff(saveat,tstops)) || dense),
              adaptive=isadaptive(alg),gamma=9//10,
              abstol=1e-2,reltol=1e-2,
              qmax=qmax_default(alg),qmin=qmin_default(alg),
              qoldinit=1//10^4, fullnormalize=true,
              beta2=beta2_default(alg),
              beta1=beta1_default(alg,beta2),
              delta=1//6,maxiters = 1e9,
              dtmax=tType((prob.tspan[end]-prob.tspan[1])),
              dtmin=tType <: AbstractFloat ? tType(10)*eps(tType) : tType(1//10^(10)),
              internalnorm=ODE_DEFAULT_NORM,
              unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
              isoutofdomain = ODE_DEFAULT_ISOUTOFDOMAIN,
              advance_to_tstop = false,stop_at_next_tstop=false,
              progress_steps=1000,
              progress=false, progress_message = ODE_DEFAULT_PROG_MESSAGE,
              progress_name="SDE",
              userdata=nothing,callback=nothing,
              timeseries_errors = true, dense_errors=false,
              kwargs...)

  tspan = prob.tspan
  tdir = sign(tspan[end]-tspan[1])

  T = tType(tspan[2])
  t = tType(tspan[1])

  if tspan[2]-tspan[1]<0 || length(tspan)>2
    error("tspan must be two numbers and final time must be greater than starting time. Aborting.")
  end

  if !(typeof(alg) <: StochasticDiffEqAdaptiveAlgorithm) && dt == 0 && isempty(tstops)
      error("Fixed timestep methods require a choice of dt or choosing the tstops")
  end

  d_discontinuities_col = collect(d_discontinuities)

  if tdir>0
    tstops_internal = binary_minheap(convert(Vector{tType},append!(collect(tstops),d_discontinuities_col)))
  else
    tstops_internal = binary_maxheap(convert(Vector{tType},append!(collect(tstops),d_discontinuities_col)))
  end

  if !isempty(tstops) && tstops[end] != tspan[2]
    push!(tstops_internal,tspan[2])
  elseif isempty(tstops)
    push!(tstops_internal,tspan[2])
  end

  if top(tstops_internal) == tspan[1]
    pop!(tstops_internal)
  end
  f = prob.f
  g = prob.g
  u0 = prob.u0
  uEltype = eltype(u0)

  (uType<:Array || uType <: Number) ? u = copy(u0) : u = deepcopy(u0)

  ks = Vector{uType}(0)

  order = alg_order(alg)

  if typeof(alg) <: SRI || typeof(alg) <: SRA
    order = alg.tableau.order
  end

  if dt == 0.0
    dt = sde_determine_initdt(u0,float(tspan[1]),tdir,dtmax,abstol,reltol,internalnorm,prob,order)
  end

  if sign(dt)!=tdir && dt!=tType(0)
    error("dt has the wrong sign. Exiting")
  end

  if typeof(u) <: AbstractArray
    rate_prototype = similar(u/zero(t),indices(u)) # rate doesn't need type info
  else
    rate_prototype = u/zero(t)
  end
  rateType = typeof(rate_prototype) ## Can be different if united

  saveat_vec =  convert(Vector{tType},collect(saveat))
  if !isempty(saveat_vec) && saveat_vec[end] == tspan[2]
    pop!(saveat_vec)
  end

  if tdir>0
    saveat_internal = binary_minheap(saveat_vec)
  else
    saveat_internal = binary_maxheap(saveat_vec)
  end

  if !isempty(saveat_internal) && top(saveat_internal) == tspan[1]
    pop!(saveat_internal)
  end

  d_discontinuities_vec =  convert(Vector{tType},d_discontinuities_col)

  if tdir>0
    d_discontinuities_internal = binary_minheap(d_discontinuities_vec)
  else
    d_discontinuities_internal = binary_maxheap(d_discontinuities_vec)
  end

  callbacks_internal = CallbackSet(callback)

  uEltypeNoUnits = typeof(recursive_one(u))
  tTypeNoUnits   = typeof(recursive_one(t))

  ### Algorithm-specific defaults ###
  ksEltype = Vector{rateType}

  # Have to convert incase passed in wrong.
  timeseries = convert(Vector{uType},timeseries_init)
  ts = convert(Vector{tType},ts_init)
  ks = convert(Vector{ksEltype},ks_init)
  alg_choice = Int[]

  copyat_or_push!(ts,1,t)
  copyat_or_push!(timeseries,1,u)
  copyat_or_push!(ks,1,[rate_prototype])

  uEltype = eltype(u)

  opts = SDEOptions(Int(maxiters),timeseries_steps,save_timeseries,adaptive,uEltype(uEltype(1)*abstol),
    uEltypeNoUnits(reltol),tTypeNoUnits(gamma),tTypeNoUnits(qmax),tTypeNoUnits(qmin),
    dtmax,dtmin,internalnorm,
    tstops_internal,saveat_internal,d_discontinuities_internal,
    userdata,
    progress,progress_steps,
    progress_name,progress_message,
    timeseries_errors,dense_errors,
    tTypeNoUnits(beta1),tTypeNoUnits(beta2),uEltypeNoUnits(delta),tTypeNoUnits(qoldinit),
    dense,save_noise,
    callbacks_internal,isoutofdomain,unstable_check,calck,advance_to_tstop,stop_at_next_tstop)

  progress ? (prog = Juno.ProgressBar(name=progress_name)) : prog = nothing

  notsaveat_idxs = Int[1]

  k = ksEltype[]

  if uType <: Array
    uprev = copy(u)
  else
    uprev = deepcopy(u)
  end


  if !(uType <: AbstractArray)
    rands = ChunkedArray(prob.noise.noise_func,ChunkedArrays.BUFFER_SIZE_DEFAULT,typeof(u/u))
    randType = typeof(u/u) # Strip units and type info
  else
    rand_prototype = similar(map((x)->x/x,u),indices(u))
    rands = ChunkedArray(prob.noise.noise_func,rand_prototype) # Strip units
    randType = typeof(rand_prototype) # Strip units and type info
  end


  Ws = Vector{randType}(0)
  if !(uType <: AbstractArray)
    W = 0.0
    Z = 0.0
    push!(Ws,W)
  else
    W = zeros(rand_prototype)
    Z = zeros(rand_prototype)
    push!(Ws,copy(W))
  end


  sqdt = sqrt(dt)
  iter = 0
  #EEst = 0
  q11 = tTypeNoUnits(1)
  ΔW = sqdt*next(rands) # Take one first
  ΔZ = sqdt*next(rands) # Take one first

  rateType = typeof(u/t) ## Can be different if united

  cache = alg_cache(alg,u,rate_prototype,ΔW,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,Val{isinplace})

  sol = build_solution(prob,alg,ts,timeseries,W=Ws,
                calculate_error = false)

  integrator =    SDEIntegrator{typeof(alg),uType,uEltype,ndims(u),ndims(u)+1,
                  tType,tTypeNoUnits,
                  uEltypeNoUnits,randType,typeof(ΔW),rateType,typeof(sol),typeof(cache),
                  typeof(prog),
                  F,F2,typeof(opts)}(f,g,uprev,t,u,tType(dt),T,alg,sol,cache,
                  rands,sqdt,W,Z,ΔW,ΔZ,opts,iter,prog,
                  tTypeNoUnits(qoldinit),q11)

  sde_solve(integrator)

  if typeof(integrator.sol.prob) <: AbstractSDETestProblem
    calculate_solution_errors!(sol;timeseries_errors=integrator.opts.timeseries_errors,dense_errors=integrator.opts.dense_errors)
  end
  sol
end
