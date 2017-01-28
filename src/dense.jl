function sde_interpolant(Θ,integrator)
  (1-Θ).*integrator.uprev .+ Θ.*integrator.u
end

function sde_interpolant(Θ,u0,u1)
  (1-Θ).*u0 .+ Θ.*u1
end

function current_interpolant(t::Number,integrator)
  Θ = (t-integrator.tprev)/integrator.dt
  sde_interpolant(Θ,integrator)
end

function current_interpolant(t::AbstractArray,integrator)
  Θ = (t.-integrator.tprev)./integrator.dt
  [sde_interpolant(ϕ,integrator) for ϕ in Θ]
end

"""
sde_interpolation(tvals,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function sde_interpolation(tvals,id)
  @unpack ts,timeseries = id
  tdir = sign(ts[end]-ts[1])
  idx = sortperm(tvals)
  i = 2 # Start the search thinking it's between ts[1] and ts[2]
  vals = Vector{eltype(timeseries)}(length(tvals))
  for j in idx
    t = tvals[j]
    i = findfirst((x)->tdir*x>=tdir*t,@view ts[i:end])+i-1 # It's in the interval ts[i-1] to ts[i]
    if ts[i] == t
      vals[j] = timeseries[i]
    elseif ts[i-1] == t # Can happen if it's the first value!
      vals[j] = timeseries[i-1]
    else
      dt = ts[i] - ts[i-1]
      Θ = (t-ts[i-1])/dt
      vals[j] = sde_interpolant(Θ,timeseries[i-1],timeseries[i])
    end
  end
  vals
end

"""
sde_interpolation(tval::Number,ts,timeseries,ks)

Get the value at tval where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function sde_interpolation(tval::Number,id)
  @unpack ts,timeseries = id
  tdir = sign(ts[end]-ts[1])
  i = findfirst((x)->tdir*x>=tdir*tval,ts) # It's in the interval ts[i-1] to ts[i]
  if ts[i] == tval
    val = timeseries[i]
  elseif ts[i-1] == tval # Can happen if it's the first value!
    push!(vals,timeseries[i-1])
  else
    dt = ts[i] - ts[i-1]
    Θ = (tval-ts[i-1])/dt
    val = sde_interpolant(Θ,timeseries[i-1],timeseries[i])
  end
  val
end
