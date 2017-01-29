@inline function change_t_via_interpolation!{T}(integrator,t,modify_save_endpoint::Type{Val{T}}=Val{false})
  # Can get rid of an allocation here with a function
  # get_tmp_arr(integrator.cache) which gives a pointer to some
  # cache array which can be modified.
  if t < integrator.tprev
    error("Current interpolant only works between tprev and t")
  elseif t != integrator.t
    integrator.dtnew = integrator.t - t
    perform_rswm_rejection!(integrator)
    new_u = integrator(t)
    if typeof(integrator.u) <: AbstractArray
      recursivecopy!(integrator.u,new_u)
    else
      integrator.u = new_u
    end
    # reeval_internals_due_to_modification!(integrator) # Not necessary for linear interp
    if T
      solution_endpoint_match_cur_integrator!(integrator)
    end
  end
end

(integrator::SDEIntegrator)(t) = current_interpolant(t,integrator)
