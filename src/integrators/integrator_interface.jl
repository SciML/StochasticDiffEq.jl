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

user_cache(integrator::SDEIntegrator) = (integrator.cache.u,integrator.cache.uprev,integrator.cache.tmp) 
u_cache(integrator::SDEIntegrator) = u_cache(integrator.cache)
du_cache(integrator::SDEIntegrator)= du_cache(integrator.cache)
full_cache(integrator::SDEIntegrator) = chain(u_cache(integrator),du_cache(integrator.cache))
