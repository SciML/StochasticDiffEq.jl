## Extrapolations are currently just constant

function sde_extrapolant!(out,Θ,integrator::DEIntegrator,idxs,deriv::Type)
  if idxs === nothing
    recursivecopy!(out,integrator.u)
  else
    out[idxs] .= @view integrator.u[idxs]
  end
end

function sde_extrapolant(Θ,integrator::DEIntegrator,idxs,deriv::Type)
  if idxs === nothing
    return integrator.u
  else
    return integrator.u[idxs]
  end
end

function sde_interpolant!(out,Θ,integrator::DEIntegrator,idxs,deriv::Type)
  sde_interpolant!(out,Θ,integrator.dt,integrator.uprev,integrator.u,idxs,deriv)
end

function sde_interpolant(Θ,integrator::DEIntegrator,idxs,deriv::Type{T}) where T
  sde_interpolant(Θ,integrator.dt,integrator.uprev,integrator.u,idxs,deriv)
end

@muladd function sde_interpolant(Θ,dt,u0,u1,idxs::Nothing,deriv::Type{Val{0}})
  @. (1-Θ)*u0 + Θ*u1
end

@muladd function sde_interpolant(Θ,dt,u0,u1,idxs,deriv::Type{Val{0}})
  @. (1-Θ)*u0[idxs] + Θ*u1[idxs]
end

function sde_interpolant(Θ,dt,u0,u1,idxs::Nothing,deriv::Type{Val{1}})
  @. (u1-u0)/dt
end

function sde_interpolant(Θ,dt,u0,u1,idxs,deriv::Type{Val{1}})
  @. (u1[idxs]-u0[idxs])/dt
end

@muladd function sde_interpolant!(out,Θ,dt,u0,u1,idxs,deriv::Type{Val{0}})
  Θm1 = (1-Θ)
  if idxs === nothing
    @. out = Θm1*u0 + Θ*u1
  else
    @views @. out = Θm1*u0[idxs] + Θ*u1[idxs]
  end
end

function sde_interpolant!(out,Θ,dt,u0,u1,idxs,deriv::Type{Val{1}})
  if idxs === nothing
    @. out = (u1-u0)/dt
  else
    @views @. out = (u1[idxs]-u0[idxs])/dt
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

@inline function current_extrapolant(t::Number,integrator::DEIntegrator,idxs=nothing,deriv=Val{0})
  Θ = (t-integrator.tprev)/(integrator.t-integrator.tprev)
  sde_extrapolant(Θ,integrator,idxs,deriv)
end

@inline function current_extrapolant!(val,t::Number,integrator::DEIntegrator,idxs=nothing,deriv=Val{0})
  Θ = (t-integrator.tprev)/(integrator.t-integrator.tprev)
  sde_extrapolant!(val,Θ,integrator,idxs,deriv)
end

@inline function current_extrapolant(t::AbstractArray,integrator::DEIntegrator,idxs=nothing,deriv=Val{0})
  Θ = (t.-integrator.tprev)./(integrator.t-integrator.tprev)
  [sde_extrapolant(ϕ,integrator,idxs,deriv) for ϕ in Θ]
end

@inline function current_extrapolant!(val,t,integrator::DEIntegrator,idxs=nothingderiv=Val{0})
  Θ = (t.-integrator.tprev)./(integrator.t-integrator.tprev)
  [sde_extrapolant!(val,ϕ,integrator,idxs,deriv) for ϕ in Θ]
end

"""
sde_interpolation(tvals,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
@inline function sde_interpolation(tvals,id,idxs,deriv,p,continuity::Symbol=:left)
  @unpack ts,timeseries = id
  tdir = sign(ts[end]-ts[1])
  idx = sortperm(tvals,rev=tdir<0)
  tdir*tvals[idx[end]] > tdir*ts[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Either solve on a longer timespan or use the local extrapolation from the integrator interface.")
  tdir*tvals[idx[1]] < tdir*ts[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Either start solving earlier or use the local extrapolation from the integrator interface.")
  i = 2 # Start the search thinking it's between ts[1] and ts[2]
  if typeof(idxs) <: Number
    vals = Vector{eltype(first(timeseries))}(undef,length(tvals))
  elseif typeof(idxs) <: AbstractArray
     vals = Vector{Array{eltype(first(timeseries)),ndims(idxs)}}(undef,length(tvals))
  else
    vals = Vector{eltype(timeseries)}(undef,length(tvals))
  end
  @inbounds for j in idx
    t = tvals[j]
    i = searchsortedfirst(@view(ts[i:end]),t,rev=tdir<0)+i-1 # It's in the interval ts[i-1] to ts[i]
    if ts[i] == t
      lasti = lastindex(ts)
      k = continuity == :right && i+1 <= lasti && ts[i+1] == t ? i+1 : i
      if idxs === nothing
        vals[j] = timeseries[k]
      else
        vals[j] = timeseries[k][idxs]
      end
    elseif ts[i-1] == t # Can happen if it's the first value!
      if idxs === nothing
        vals[j] = timeseries[i-1]
      else
        vals[j] = timeseries[i-1][idxs]
      end
    else
      dt = ts[i] - ts[i-1]
      Θ = (t-ts[i-1])/dt
      vals[j] = sde_interpolant(Θ,dt,timeseries[i-1],timeseries[i],idxs,deriv)
    end
  end
  DiffEqArray(vals,tvals)
end

"""
sde_interpolation(tval::Number,ts,timeseries,ks)

Get the value at tval where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
@inline function sde_interpolation(tval::Number,id,idxs,deriv,p,continuity::Symbol=:left)
  @unpack ts,timeseries = id
  tdir = sign(ts[end]-ts[1])
  tdir*tval > tdir*ts[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Either solve on a longer timespan or use the local extrapolation from the integrator interface.")
  tdir*tval < tdir*ts[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Either start solving earlier or use the local extrapolation from the integrator interface.")
  @inbounds i = searchsortedfirst(ts,tval,rev=tdir<0) # It's in the interval ts[i-1] to ts[i]
  @inbounds if ts[i] == tval
    lasti = lastindex(ts)
    k = continuity == :right && i+1 <= lasti && ts[i+1] == tval ? i+1 : i
    if idxs === nothing
      val = timeseries[k]
    else
      val = timeseries[k][idxs]
    end
  elseif ts[i-1] == tval # Can happen if it's the first value!
    if idxs === nothing
      val = timeseries[i-1]
    else
      val = timeseries[i-1][idxs]
    end
  else
    dt = ts[i] - ts[i-1]
    Θ = (tval-ts[i-1])/dt
    val = sde_interpolant(Θ,dt,timeseries[i-1],timeseries[i],idxs,deriv)
  end
  val
end

@inline function sde_interpolation!(out,tval::Number,id,idxs,deriv,p,continuity::Symbol=:left)
  @unpack ts,timeseries = id
  tdir = sign(ts[end]-ts[1])
  tdir*tval > tdir*ts[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Either solve on a longer timespan or use the local extrapolation from the integrator interface.")
  tdir*tval < tdir*ts[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Either start solving earlier or use the local extrapolation from the integrator interface.")
  @inbounds i = searchsortedfirst(ts,tval,rev=tdir<0) # It's in the interval ts[i-1] to ts[i]
  @inbounds if ts[i] == tval
    lasti = lastindex(ts)
    k = continuity == :right && i+1 <= lasti && ts[i+1] == tval ? i+1 : i
    if idxs === nothing
      copyto!(out,timeseries[k])
    else
      copyto!(out,timeseries[k][idxs])
    end
  elseif ts[i-1] == tval # Can happen if it's the first value!
    if idxs === nothing
      copyto!(out,timeseries[i-1])
    else
      copyto!(out,timeseries[i-1][idxs])
    end
  else
    dt = ts[i] - ts[i-1]
    Θ = (tval-ts[i-1])/dt
    sde_interpolant!(out,Θ,dt,timeseries[i-1],timeseries[i],idxs,deriv)
  end
end

@inline function sde_interpolation!(vals,tvals,id,idxs,deriv,p,continuity::Symbol=:left)
  @unpack ts,timeseries = id
  tdir = sign(ts[end]-ts[1])
  idx = sortperm(tvals,rev=tdir<0)
  tdir*tvals[idx[end]] > tdir*ts[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Either solve on a longer timespan or use the local extrapolation from the integrator interface.")
  tdir*tvals[idx[1]] < tdir*ts[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Either start solving earlier or use the local extrapolation from the integrator interface.")
  i = 2 # Start the search thinking it's between ts[1] and ts[2]
  @inbounds for j in idx
    t = tvals[j]
    i = searchsortedfirst(@view(ts[i:end]),t,rev=tdir<0)+i-1 # It's in the interval ts[i-1] to ts[i]
    if ts[i] == t
      lasti = lastindex(ts)
      k = continuity == :right && i+1 <= lasti && ts[i+1] == t ? i+1 : i
      if idxs === nothing
        vals[j] = timeseries[k]
      else
        vals[j] = timeseries[k][idxs]
      end
    elseif ts[i-1] == t # Can happen if it's the first value!
      if idxs === nothing
        vals[j] = timeseries[i-1]
      else
        vals[j] = timeseries[i-1][idxs]
      end
    else
      dt = ts[i] - ts[i-1]
      Θ = (t-ts[i-1])/dt
      if eltype(timeseries) <: AbstractArray
        sde_interpolant!(vals[j],Θ,dt,timeseries[i-1],timeseries[i],idxs,deriv)
      else
        vals[j] = sde_interpolant(Θ,dt,timeseries[i-1],timeseries[i],idxs,deriv)
      end
    end
  end
end
