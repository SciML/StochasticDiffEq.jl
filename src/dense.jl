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
function sde_interpolation(tvals,sol)
  @unpack ts,timeseries = sol
  tdir = sign(ts[end]-ts[1])
  idx = sortperm(tvals)
  i = 2 # Start the search thinking it's between ts[1] and ts[2]
  vals = Vector{eltype(timeseries)}(length(tvals))
  for j in idx
    t = tvals[j]
    i = findfirst((x)->tdir*x>=tdir*t,ts[notsaveat_idxs[i:end]])+i-1 # It's in the interval ts[i-1] to ts[i]
    if ts[notsaveat_idxs[i]] == t
      vals[j] = timeseries[notsaveat_idxs[i]]
    elseif ts[notsaveat_idxs[i-1]] == t # Can happen if it's the first value!
      vals[j] = timeseries[notsaveat_idxs[i-1]]
    else
      dt = ts[notsaveat_idxs[i]] - ts[notsaveat_idxs[i-1]]
      Θ = (t-ts[notsaveat_idxs[i-1]])/dt
      vals[j] = sde_interpolant(Θ,timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]])
    end
  end
  vals
end

"""
sde_interpolation(tval::Number,ts,timeseries,ks)

Get the value at tval where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function ode_interpolation(tval::Number,sol)
  @unpack ts,timeseries = sol
  tdir = sign(ts[end]-ts[1])
  i = findfirst((x)->tdir*x>=tdir*tval,@view ts[notsaveat_idxs]) # It's in the interval ts[i-1] to ts[i]
  if ts[notsaveat_idxs[i]] == tval
    val = timeseries[notsaveat_idxs[i]]
  elseif ts[notsaveat_idxs[i-1]] == tval # Can happen if it's the first value!
    push!(vals,timeseries[notsaveat_idxs[i-1]])
  else
    dt = ts[notsaveat_idxs[i]] - ts[notsaveat_idxs[i-1]]
    Θ = (tval-ts[notsaveat_idxs[i-1]])/dt
    val = sde_interpolant(Θ,timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]])
  end
  val
end
