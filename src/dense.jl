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
  @.. (1-Θ)*u0 + Θ*u1
end

@muladd function sde_interpolant(Θ,dt,u0,u1,idxs,deriv::Type{Val{0}})
  @.. (1-Θ)*u0[idxs] + Θ*u1[idxs]
end

function sde_interpolant(Θ,dt,u0,u1,idxs::Nothing,deriv::Type{Val{1}})
  @.. (u1-u0)/dt
end

function sde_interpolant(Θ,dt,u0,u1,idxs,deriv::Type{Val{1}})
  @.. (u1[idxs]-u0[idxs])/dt
end

@muladd function sde_interpolant!(out,Θ,dt,u0,u1,idxs,deriv::Type{Val{0}})
  Θm1 = (1-Θ)
  if idxs === nothing
    @.. out = Θm1*u0 + Θ*u1
  else
    @views @.. out = Θm1*u0[idxs] + Θ*u1[idxs]
  end
end

function sde_interpolant!(out,Θ,dt,u0,u1,idxs,deriv::Type{Val{1}})
  if idxs === nothing
    @.. out = (u1-u0)/dt
  else
    @views @.. out = (u1[idxs]-u0[idxs])/dt
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
function sde_interpolation(tvals,id,idxs,deriv,p,continuity::Symbol=:left)
  @unpack ts,timeseries = id
  @inbounds tdir = sign(ts[end]-ts[1])
  idx = sortperm(tvals,rev=tdir<0)

  # start the search thinking it's ts[1]-ts[2]
  i₋ = 1
  i₊ = 2
  vals = map(idx) do j
    @inbounds begin
      t = tvals[j]

      if continuity === :left
        # we have i₋ = i₊ = 1 if t = ts[1], i₊ = i₋ + 1 = lastindex(ts) if t > ts[end],
        # and otherwise i₋ and i₊ satisfy ts[i₋] < t ≤ ts[i₊]
        i₊ = min(lastindex(ts), OrdinaryDiffEq._searchsortedfirst(ts,t,i₊,tdir > 0))
        i₋ = i₊ > 1 ? i₊ - 1 : i₊
      else
        # we have i₋ = i₊ - 1 = 1 if t < ts[1], i₊ = i₋ = lastindex(ts) if t = ts[end],
        # and otherwise i₋ and i₊ satisfy ts[i₋] ≤ t < ts[i₊]
        i₋ = max(1, OrdinaryDiffEq._searchsortedlast(ts,t,i₋,tdir > 0))
        i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
      end

      dt = ts[i₊] - ts[i₋]
      Θ = iszero(dt) ? oneunit(t) / oneunit(dt) : (t-ts[i₋]) / dt

      sde_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],idxs,deriv)
    end
  end

  DiffEqArray(vals, tvals)
end

"""
sde_interpolation(tval::Number,ts,timeseries,ks)

Get the value at tval where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function sde_interpolation(tval::Number,id,idxs,deriv,p,continuity::Symbol=:left)
  @unpack ts,timeseries = id
  @inbounds tdir = sign(ts[end]-ts[1])

  if continuity === :left
    # we have i₋ = i₊ = 1 if tval = ts[1], i₊ = i₋ + 1 = lastindex(ts) if tval > ts[end],
    # and otherwise i₋ and i₊ satisfy ts[i₋] < tval ≤ ts[i₊]
    i₊ = min(lastindex(ts), OrdinaryDiffEq._searchsortedfirst(ts,tval,2,tdir > 0))
    i₋ = i₊ > 1 ? i₊ - 1 : i₊
  else
    # we have i₋ = i₊ - 1 = 1 if tval < ts[1], i₊ = i₋ = lastindex(ts) if tval = ts[end],
    # and otherwise i₋ and i₊ satisfy ts[i₋] ≤ tval < ts[i₊]
    i₋ = max(1, OrdinaryDiffEq._searchsortedlast(ts,tval,1,tdir > 0))
    i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
  end

  @inbounds begin
    dt = ts[i₊] - ts[i₋]
    Θ = iszero(dt) ? oneunit(tval) / oneunit(dt) : (tval-ts[i₋]) / dt
    val = sde_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],idxs,deriv)
  end

  val
end

function sde_interpolation!(out,tval::Number,id,idxs,deriv,p,continuity::Symbol=:left)
  @unpack ts,timeseries = id
  @inbounds tdir = sign(ts[end]-ts[1])

  if continuity === :left
    # we have i₋ = i₊ = 1 if tval = ts[1], i₊ = i₋ + 1 = lastindex(ts) if tval > ts[end],
    # and otherwise i₋ and i₊ satisfy ts[i₋] < tval ≤ ts[i₊]
    i₊ = min(lastindex(ts), OrdinaryDiffEq._searchsortedfirst(ts,tval,2,tdir > 0))
    i₋ = i₊ > 1 ? i₊ - 1 : i₊
  else
    # we have i₋ = i₊ - 1 = 1 if tval < ts[1], i₊ = i₋ = lastindex(ts) if tval = ts[end],
    # and otherwise i₋ and i₊ satisfy ts[i₋] ≤ tval < ts[i₊]
    i₋ = max(1, OrdinaryDiffEq._searchsortedlast(ts,tval,1,tdir > 0))
    i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
  end

  @inbounds begin
    dt = ts[i₊] - ts[i₋]
    Θ = iszero(dt) ? oneunit(tval) / oneunit(dt) : (tval-ts[i₋]) / dt
    sde_interpolant!(out,Θ,dt,timeseries[i₋],timeseries[i₊],idxs,deriv)
  end

  out
end

function sde_interpolation!(vals,tvals,id,idxs,deriv,p,continuity::Symbol=:left)
  @unpack ts,timeseries = id
  @inbounds tdir = sign(ts[end]-ts[1])
  idx = sortperm(tvals,rev=tdir<0)

  # start the search thinking it's in ts[1]-ts[2]
  i₋ = 1
  i₊ = 2
  @inbounds for j in idx
    t = tvals[j]

    if continuity === :left
      # we have i₋ = i₊ = 1 if t = ts[1], i₊ = i₋ + 1 = lastindex(ts) if t > ts[end],
      # and otherwise i₋ and i₊ satisfy ts[i₋] < t ≤ ts[i₊]
      i₊ = min(lastindex(ts), OrdinaryDiffEq._searchsortedfirst(ts,t,i₊,tdir > 0))
      i₋ = i₊ > 1 ? i₊ - 1 : i₊
    else
      # we have i₋ = i₊ - 1 = 1 if t < ts[1], i₊ = i₋ = lastindex(ts) if t = ts[end],
      # and otherwise i₋ and i₊ satisfy ts[i₋] ≤ t < ts[i₊]
      i₋ = max(1, OrdinaryDiffEq._searchsortedlast(ts,t,i₋,tdir > 0))
      i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
    end

    dt = ts[i₊] - ts[i₋]
    Θ = iszero(dt) ? oneunit(t) / oneunit(dt) : (t-ts[i₋]) / dt

    if eltype(timeseries) <: AbstractArray
      sde_interpolant!(vals[j],Θ,dt,timeseries[i₋],timeseries[i₊],idxs,deriv)
    else
      vals[j] = sde_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],idxs,deriv)
    end
  end

  vals
end
