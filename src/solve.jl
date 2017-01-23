@inline ODE_DEFAULT_NORM(u) = sqrt(sum(abs2,u) / length(u))
@inline ODE_DEFAULT_PROG_MESSAGE(dt,t,u) = "dt="*string(dt)*"\nt="*string(t)*"\nmax u="*string(maximum(abs.(u)))
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,t,u) = any(isnan,u)

function solve{uType,tType,isinplace,NoiseClass,F,F2,F3,algType<:AbstractSDEAlgorithm,recompile_flag}(
              prob::AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3},
              alg::algType,timeseries=[],ts=[],ks=[],recompile::Type{Val{recompile_flag}}=Val{true};
              dt = tType(0),save_timeseries::Bool = true,
              timeseries_steps::Int = 1,
              dense = false,
              saveat = tType[],tstops = tType[],d_discontinuities= tType[],
              calck = (!isempty(setdiff(saveat,tstops)) || dense),
              adaptive=isadaptive(alg),γ=9//10,
              abstol=1e-2,reltol=1e-2,
              qmax=qmax_default(alg),qmin=qmin_default(alg),
              qoldinit=1//10^4, fullnormalize=true,
              beta2=beta2_default(alg),
              beta1=beta1_default(alg,beta2),
              δ=1/6,maxiters = 1e9,
              dtmax=tType((prob.tspan[end]-prob.tspan[1])),
              dtmin=tType <: AbstractFloat ? tType(10)*eps(tType) : tType(1//10^(10)),
              internalnorm=ODE_DEFAULT_NORM,
              unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
              advance_to_tstop = false,stop_at_next_tstop=false,
              discard_length=1e-15,adaptivealg::Symbol=:RSwM3,
              progress_steps=1000,
              progress=false, progress_message = ODE_DEFAULT_PROG_MESSAGE,
              progress_name="SDE",
              userdata=nothing,callback=nothing,
              timeseries_errors = true, dense_errors=false,
              kwargs...)

  @unpack u0,noise,tspan = prob

  tspan = prob.tspan
  tdir = sign(tspan[end]-tspan[1])

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

  T = tType(tspan[2])
  t = tType(tspan[1])

  timeseries = Vector{uType}(0)
  push!(timeseries,u0)
  ts = Vector{tType}(0)
  push!(ts,t)

  uEltype = eltype(u)

  if !(uType <: AbstractArray)
    rands = ChunkedArray(noise.noise_func)
    randType = typeof(u/u) # Strip units and type info
  else
    rand_prototype = similar(map((x)->x/x,u),indices(u))
    rands = ChunkedArray(noise.noise_func,rand_prototype) # Strip units
    randType = typeof(rand_prototype) # Strip units and type info
  end

  uEltypeNoUnits = typeof(recursive_one(u))
  tTypeNoUnits   = typeof(recursive_one(t))


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
  maxstacksize = 0
  #EEst = 0
  q11 = tTypeNoUnits(1)

  rateType = typeof(u/t) ## Can be different if united

  #@code_warntype sde_solve(SDEIntegrator{typeof(alg),typeof(u),eltype(u),ndims(u),ndims(u)+1,typeof(dt),typeof(tableau)}(f,g,u,t,dt,T,Int(maxiters),timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,dtmax,dtmin,internalnorm,discard_length,progress,atomloaded,progress_steps,rands,sqdt,W,Z,tableau))

  u,t,W,timeseries,ts,Ws,maxstacksize,maxstacksize2 = sde_solve(
  SDEIntegrator{typeof(alg),uType,uEltype,ndims(u),ndims(u)+1,tType,tTypeNoUnits,
                uEltypeNoUnits,randType,rateType,typeof(internalnorm),typeof(progress_message),
                typeof(unstable_check),F,F2}(f,g,u,t,dt,T,alg,Int(maxiters),timeseries,Ws,
                ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,tTypeNoUnits(γ),
                abstol,reltol,tTypeNoUnits(qmax),dtmax,dtmin,internalnorm,discard_length,
                progress,progress_name,progress_steps,progress_message,
                unstable_check,rands,sqdt,W,Z,
                tTypeNoUnits(beta1),tTypeNoUnits(beta2),tTypeNoUnits(qoldinit),tTypeNoUnits(qmin),q11,tTypeNoUnits(qoldinit)))

  build_solution(prob,alg,ts,timeseries,W=Ws,
                  timeseries_errors = timeseries_errors,
                  maxstacksize = maxstacksize)

end
