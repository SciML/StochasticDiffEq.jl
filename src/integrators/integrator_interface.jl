@inline function change_t_via_interpolation!{T}(integrator,t,modify_save_endpoint::Type{Val{T}}=Val{false})
  # Can get rid of an allocation here with a function
  # get_tmp_arr(integrator.cache) which gives a pointer to some
  # cache array which can be modified.
  if t < integrator.tprev
    error("Current interpolant only works between tprev and t")
  elseif t != integrator.t
    if typeof(integrator.u) <: AbstractArray
      integrator(integrator.u,t)
    else
      integrator.u = integrator(t)
    end
    integrator.dtnew = integrator.t - t
    perform_rswm_rejection!(integrator) #this only changes dt and noise, so no interpolation problems
    integrator.t = t
    # reeval_internals_due_to_modification!(integrator) # Not necessary for linear interp
    if T
      solution_endpoint_match_cur_integrator!(integrator)
    end
  end
end

function (integrator::SDEIntegrator)(t,deriv::Type=Val{0};idxs=size(integrator.uprev))
  current_interpolant(t,integrator,idxs,deriv)
end

(integrator::SDEIntegrator)(val::AbstractArray,t::Union{Number,AbstractArray},deriv::Type=Val{0};idxs=eachindex(integrator.uprev)) = current_interpolant!(val,t,integrator,idxs,deriv)


user_cache(integrator::SDEIntegrator) = user_cache(integrator.cache)
u_cache(integrator::SDEIntegrator) = u_cache(integrator.cache)
du_cache(integrator::SDEIntegrator)= du_cache(integrator.cache)
user_cache(c::StochasticDiffEqCache) = (c.u,c.uprev,c.tmp)
full_cache(integrator::SDEIntegrator) = chain(user_cache(integrator),u_cache(integrator),du_cache(integrator.cache))
default_non_user_cache(integrator::SDEIntegrator) = chain(u_cache(integrator),du_cache(integrator.cache))
@inline add_tstop!(integrator::SDEIntegrator,t) = push!(integrator.opts.tstops,t)

resize_non_user_cache!(integrator::SDEIntegrator,i::Int) = resize_non_user_cache!(integrator,integrator.cache,i)
resize!(integrator::SDEIntegrator,i::Int) = resize!(integrator,integrator.cache,i)

function resize!(integrator::SDEIntegrator,cache,i)
  resize_non_user_cache!(integrator,cache,i)
  for c in user_cache(integrator)
    resize!(c,i)
  end
end

function resize_noise!(integrator,cache,prev_len,i)
  for c in integrator.S₁
    resize!(c[2],i)
    resize!(c[3],i)
    if i > prev_len # fill in rands
      resize_noise_caches!(integrator,c,c[1],prev_len:i)
    end
  end
  for c in integrator.S₂
    resize!(c[2],i)
    resize!(c[3],i)
    if i > prev_len # fill in rands
      resize_noise_caches!(integrator,c,c[1],prev_len:i)
    end
  end
  resize!(integrator.ΔW,i)
  resize!(integrator.ΔZ,i)
  resize!(integrator.ΔWtilde,i)
  resize!(integrator.ΔZtilde,i)
  resize!(integrator.ΔWtmp,i)
  resize!(integrator.ΔZtmp,i)
  resize!(integrator.W,i)
  resize!(integrator.Z,i)
  if i > prev_len # fill in rands
    fill!(@view(integrator.W[prev_len:i]),zero(eltype(integrator.u)))
    fill!(@view(integrator.Z[prev_len:i]),zero(eltype(integrator.u)))
  end
end

function resize_non_user_cache!(integrator::SDEIntegrator,cache,i)
  prev_len = length(integrator.u)
  resize_noise!(integrator,cache,prev_len,i)
  for c in default_non_user_cache(integrator)
    resize!(c,i)
  end
end

function deleteat!(integrator::SDEIntegrator,i::Int)
  for c in full_cache(integrator)
    deleteat!(c,i)
  end
end

function terminate!(integrator::SDEIntegrator)
  integrator.opts.tstops.valtree = typeof(integrator.opts.tstops.valtree)()
end
