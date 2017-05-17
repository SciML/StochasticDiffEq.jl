## Extrapolations are currently just constant

function sde_extrapolant!(out,Θ,integrator::DEIntegrator,idxs,deriv::Type)
  out .= integrator.u
end

function sde_extrapolant(Θ,integrator::DEIntegrator,idxs,deriv::Type)
  integrator.u
end

function sde_interpolant!(out,Θ,integrator::DEIntegrator,idxs,deriv::Type)
  sde_interpolant!(out,Θ,integrator.dt,integrator.uprev,integrator.u,idxs,deriv)
end

function sde_interpolant(Θ,integrator::DEIntegrator,idxs,deriv::Type)
  sde_interpolant(Θ,integrator.dt,integrator.uprev,integrator.u,idxs,deriv)
end

function sde_interpolant(Θ,dt,u0::Number,u1,idxs,deriv::Type{Val{0}})
  (1-Θ)*u0 + Θ*u1
end

function sde_interpolant(Θ,dt,u0::Number,u1,idxs,deriv::Type{Val{1}})
  (u1-u0)/dt
end

function sde_interpolant!(out,Θ,dt,u0,u1,idxs,deriv::Type{Val{0}})
  Θm1 = (1-Θ)
  for (j,i) in enumerate(idxs)
    out[j] = Θm1*u0[i] + Θ*u1[i]
  end
end

function sde_interpolant!(out,Θ,dt,u0,u1,idxs,deriv::Type{Val{1}})
  for (j,i) in enumerate(idxs)
    out[j] = (u1[i]-u0[i])/dt
  end
end

function sde_interpolant(Θ,dt,u0::AbstractArray,u1,idxs,deriv::Type)
  if typeof(idxs) <: Tuple
    out = similar(u0,idxs)
    idxs_internal=eachindex(u0)
  else
    out = similar(u0,indices(idxs))
    idxs_internal=idxs
  end
  sde_interpolant!(out,Θ,dt,u0,u1,idxs_internal,deriv)
  if typeof(idxs) <: Number
    return first(out)
  else
    return out
  end
end

@inline function current_interpolant(t::Number,integrator::DEIntegrator,idxs,deriv)
  Θ = (t-integrator.tprev)/integrator.dt
  sde_interpolant(Θ,integrator,idxs,deriv)
end

@inline function current_interpolant(t,integrator::DEIntegrator,idxs,deriv)
  Θ = (t.-integrator.tprev)./integrator.dt
  [sde_interpolant(ϕ,integrator,idxs,deriv) for ϕ in Θ]
end

@inline function current_interpolant!(val,t::Number,integrator::DEIntegrator,idxs,deriv)
  Θ = (t-integrator.tprev)/integrator.dt
  sde_interpolant!(val,Θ,integrator,idxs,deriv)
end

@inline function current_interpolant!(val,t,integrator::DEIntegrator,idxs,deriv)
  Θ = (t.-integrator.tprev)./integrator.dt
  [sde_interpolant!(val,ϕ,integrator,idxs,deriv) for ϕ in Θ]
end

@inline function current_extrapolant(t::Number,integrator::DEIntegrator,idxs=size(integrator.uprev),deriv=Val{0})
  Θ = (t-integrator.tprev)/(integrator.t-integrator.tprev)
  sde_extrapolant(Θ,integrator,idxs,deriv)
end

@inline function current_extrapolant!(val,t::Number,integrator::DEIntegrator,idxs=size(integrator.uprev),deriv=Val{0})
  Θ = (t-integrator.tprev)/(integrator.t-integrator.tprev)
  sde_extrapolant!(val,Θ,integrator,idxs,deriv)
end

@inline function current_extrapolant(t::AbstractArray,integrator::DEIntegrator,idxs=size(integrator.uprev),deriv=Val{0})
  Θ = (t.-integrator.tprev)./(integrator.t-integrator.tprev)
  [sde_extrapolant(ϕ,integrator,idxs,deriv) for ϕ in Θ]
end

@inline function current_extrapolant!(val,t,integrator::DEIntegrator,idxs=size(integrator.uprev),deriv=Val{0})
  Θ = (t.-integrator.tprev)./(integrator.t-integrator.tprev)
  [sde_extrapolant!(val,ϕ,integrator,idxs,deriv) for ϕ in Θ]
end

"""
sde_interpolation(tvals,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
@inline function sde_interpolation(tvals,id,idxs,deriv)
  @unpack ts,timeseries = id
  tdir = sign(ts[end]-ts[1])
  idx = sortperm(tvals,rev=tdir<0)
  tvals[idx[end]] > ts[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Either solve on a longer timespan or use the local extrapolation from the integrator interface.")
  tvals[idx[1]] < ts[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Either start solving earlier or use the local extrapolation from the integrator interface.")
  i = 2 # Start the search thinking it's between ts[1] and ts[2]
  if idxs == nothing
    if (eltype(timeseries) <: AbstractArray) && !(eltype(timeseries) <: Union{StaticArray,Array})
      vals = Vector{Vector{eltype(first(timeseries))}}(length(tvals))
    else
      vals = Vector{eltype(timeseries)}(length(tvals))
    end
  elseif typeof(idxs) <: Number
    vals = Vector{eltype(first(timeseries))}(length(tvals))
  else
    vals = Vector{Vector{eltype(first(timeseries))}}(length(tvals))
  end
  @inbounds for j in idx
    t = tvals[j]
    i = searchsortedfirst(@view(ts[i:end]),t,rev=tdir<0)+i-1 # It's in the interval ts[i-1] to ts[i]
    if ts[i] == t
      if idxs == nothing
        vals[j] = timeseries[i]
      else
        vals[j] = timeseries[i][idxs]
      end
    elseif ts[i-1] == t # Can happen if it's the first value!
      if idxs == nothing
        vals[j] = timeseries[i-1]
      else
        vals[j] = timeseries[i-1][idxs]
      end
    else
      dt = ts[i] - ts[i-1]
      Θ = (t-ts[i-1])/dt
      if idxs == nothing && eltype(timeseries) <: AbstractArray
        idxs_internal = size(timeseries[i-1])
      else
        idxs_internal = idxs
      end
      vals[j] = sde_interpolant(Θ,dt,timeseries[i-1],timeseries[i],idxs_internal,deriv)
    end
  end
  vals
end

"""
sde_interpolation(tval::Number,ts,timeseries,ks)

Get the value at tval where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
@inline function sde_interpolation(tval::Number,id,idxs,deriv)
  @unpack ts,timeseries = id
  tval > ts[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Either solve on a longer timespan or use the local extrapolation from the integrator interface.")
  tval < ts[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Either start solving earlier or use the local extrapolation from the integrator interface.")
  tdir = sign(ts[end]-ts[1])
  @inbounds i = searchsortedfirst(ts,tval,rev=tdir<0) # It's in the interval ts[i-1] to ts[i]
  @inbounds if ts[i] == tval
    if idxs == nothing
      val = timeseries[i]
    else
      val = timeseries[i][idxs]
    end
  elseif ts[i-1] == tval # Can happen if it's the first value!
    if idxs == nothing
      val = timeseries[i-1]
    else
      val = timeseries[i-1][idxs]
    end
  else
    dt = ts[i] - ts[i-1]
    Θ = (tval-ts[i-1])/dt
    if idxs == nothing && eltype(timeseries) <: AbstractArray
      idxs_internal = size(timeseries[i-1])
    else
      idxs_internal = idxs
    end
    val = sde_interpolant(Θ,dt,timeseries[i-1],timeseries[i],idxs_internal,deriv)
  end
  val
end

@inline function sde_interpolation!(out,tval::Number,id,idxs,deriv)
  @unpack ts,timeseries = id
  tval > ts[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Either solve on a longer timespan or use the local extrapolation from the integrator interface.")
  tval < ts[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Either start solving earlier or use the local extrapolation from the integrator interface.")
  tdir = sign(ts[end]-ts[1])
  @inbounds i = searchsortedfirst(ts,tval,rev=tdir<0) # It's in the interval ts[i-1] to ts[i]
  @inbounds if ts[i] == tval
    if idxs == nothing
      copy!(out,timeseries[i])
    else
      copy!(out,timeseries[i][idxs])
    end
  elseif ts[i-1] == tval # Can happen if it's the first value!
    if idxs == nothing
      copy!(out,timeseries[i-1])
    else
      copy!(out,timeseries[i-1][idxs])
    end
  else
    dt = ts[i] - ts[i-1]
    Θ = (tval-ts[i-1])/dt
    if idxs == nothing
      idxs_internal = eachindex(out)
    else
      idxs_internal = idxs
    end
    sde_interpolant!(out,Θ,dt,timeseries[i-1],timeseries[i],idxs_internal,deriv)
  end
end

@inline function sde_interpolation!(vals,tvals,id,idxs,deriv)
  @unpack ts,timeseries = id
  tdir = sign(ts[end]-ts[1])
  idx = sortperm(tvals,rev=tdir<0)
  tvals[idx[end]] > ts[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Either solve on a longer timespan or use the local extrapolation from the integrator interface.")
  tvals[idx[1]] < ts[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Either start solving earlier or use the local extrapolation from the integrator interface.")
  i = 2 # Start the search thinking it's between ts[1] and ts[2]
  @inbounds for j in idx
    t = tvals[j]
    i = searchsortedfirst(@view(ts[i:end]),t,rev=tdir<0)+i-1 # It's in the interval ts[i-1] to ts[i]
    if ts[i] == t
      if idxs == nothing
        vals[j] = timeseries[i]
      else
        vals[j] = timeseries[i][idxs]
      end
    elseif ts[i-1] == t # Can happen if it's the first value!
      if idxs == nothing
        vals[j] = timeseries[i-1]
      else
        vals[j] = timeseries[i-1][idxs]
      end
    else
      dt = ts[i] - ts[i-1]
      Θ = (t-ts[i-1])/dt
      if idxs == nothing && eltype(vals) <: AbstractArray
        idxs_internal = eachindex(vals[j])
      else
        idxs_internal = idxs
      end
      if eltype(timeseries) <: AbstractArray
        sde_interpolant!(vals[j],Θ,dt,timeseries[i-1],timeseries[i],idxs_internal,deriv)
      else
        vals[j] = sde_interpolant(Θ,dt,timeseries[i-1],timeseries[i],idxs_internal,deriv)
      end
    end
  end
end
