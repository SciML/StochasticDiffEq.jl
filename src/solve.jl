@inline ODE_DEFAULT_NORM(u) = sqrt(sum(abs2,u) / length(u))
@inline ODE_DEFAULT_PROG_MESSAGE(dt,t,u) = "dt="*string(dt)*"\nt="*string(t)*"\nmax u="*string(maximum(abs.(u)))
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,t,u) = any(isnan,u)

function solve{uType,tType,isinplace,NoiseClass,F,F2,F3,algType<:AbstractSDEAlgorithm}(
              prob::AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3},
              alg::algType;
              dt::Number=0.0,save_timeseries::Bool = true,
              timeseries_steps::Int = 1,adaptive=true,γ=2.0,alg_hint=nothing,
              abstol=1e-2,reltol=1e-2,qmax=1.125,δ=1/6,maxiters::Int = round(Int,1e9),
              dtmax=nothing,dtmin=nothing,internalnorm=ODE_DEFAULT_NORM,
              tstops=tType[],saveat=tType[],
              unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
              discard_length=1e-15,adaptivealg::Symbol=:RSwM3,
              progress_steps=1000,
              progress=false, progress_message = ODE_DEFAULT_PROG_MESSAGE,
              progress_name="SDE",
              timeseries_errors=true,tableau = nothing,kwargs...)

  @unpack u0,noise,tspan = prob

  tdir = tspan[2]-tspan[1]

  if tspan[2]-tspan[1]<0 || length(tspan)>2
    error("tspan must be two numbers and final time must be greater than starting time. Aborting.")
  end

  if !(typeof(alg) <: StochasticDiffEqAdaptiveAlgorithm) && dt == 0 && isempty(tstops)
      error("Fixed timestep methods require a choice of dt or choosing the tstops")
  end

  u = copy(u0)
  if !isinplace && typeof(u)<:AbstractArray
    f = (t,u,du) -> (du[:] = prob.f(t,u))
    g = (t,u,du) -> (du[:] = prob.g(t,u))
  else
    f = prob.f
    g = prob.g
  end

  if typeof(alg) <: StochasticDiffEqAdaptiveAlgorithm
    if adaptive
      dt = 1.0*dt
    end
  else
    adaptive = false
  end

  if dtmax == nothing
    dtmax = (tspan[2]-tspan[1])/2
  end
  if dtmin == nothing
    if tType <: AbstractFloat
      dtmin = tType(10)*eps(tType)
    else
      dtmin = tType(1//10^(10))
    end
  end

  if dt == 0.0
    if typeof(alg)==EM
      order = 0.5
    elseif typeof(alg)==RKMil
      order = 1.0
    elseif typeof(alg) == SRA1 || typeof(alg) == SRA
      order = 2.0
    else
      order = 1.5
    end
    dt = sde_determine_initdt(u0,float(tspan[1]),tdir,dtmax,abstol,reltol,internalnorm,f,g,order)
  end

  T = tType(tspan[2])
  t = tType(tspan[1])

  timeseries = Vector{uType}(0)
  push!(timeseries,u0)
  ts = Vector{tType}(0)
  push!(ts,t)

  #PreProcess
  if typeof(alg)== SRA && tableau == nothing
    tableau = constructSRA1()
  elseif tableau == nothing # && (typeof(alg)==:SRI)
    tableau = constructSRIW1()
  end

  uEltype = eltype(u)
  tableauType = typeof(tableau)

  if !(uType <: AbstractArray)
    rands = ChunkedArray(noise.noise_func)
  else
    rands = ChunkedArray(noise.noise_func,map((x)->x/x,u)) # Strip units
  end

  randType = typeof(map((x)->x/x,u)) # Strip units

  if uType <: AbstractArray
    uEltypeNoUnits = eltype(u./u)
  else
    uEltypeNoUnits = typeof(u./u)
  end


  Ws = Vector{randType}(0)
  if !(uType <: AbstractArray)
    W = 0.0
    Z = 0.0
    push!(Ws,W)
  else
    W = zeros(u0)
    Z = zeros(u0)
    push!(Ws,copy(W))
  end
  sqdt = sqrt(dt)
  iter = 0
  maxstacksize = 0
  #EEst = 0

  rateType = typeof(u/t) ## Can be different if united

  #@code_warntype sde_solve(SDEIntegrator{typeof(alg),typeof(u),eltype(u),ndims(u),ndims(u)+1,typeof(dt),typeof(tableau)}(f,g,u,t,dt,T,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,dtmax,dtmin,internalnorm,discard_length,progress,atomloaded,progress_steps,rands,sqdt,W,Z,tableau))

  u,t,W,timeseries,ts,Ws,maxstacksize,maxstacksize2 = sde_solve(SDEIntegrator{typeof(alg),uType,uEltype,ndims(u),ndims(u)+1,tType,tableauType,uEltypeNoUnits,randType,rateType,typeof(internalnorm),typeof(progress_message),typeof(unstable_check),typeof(f),typeof(g)}(f,g,u,t,dt,T,alg,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,dtmax,dtmin,internalnorm,discard_length,progress,progress_name,progress_steps,progress_message,unstable_check,rands,sqdt,W,Z,tableau))

  build_solution(prob,alg,ts,timeseries,W=Ws,
                  timeseries_errors = timeseries_errors,
                  maxstacksize = maxstacksize)

end
