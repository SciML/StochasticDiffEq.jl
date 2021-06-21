@inline function initialize!(integrator,cache::StochasticCompositeCache,f::F=integrator.f) where F
  cache.current = cache.choice_function(integrator)
  if cache.current == 1
    initialize!(integrator, @inbounds(cache.caches[1]))
  elseif cache.current == 2
    initialize!(integrator, @inbounds(cache.caches[2]))
  else
    initialize!(integrator, @inbounds(cache.caches[cache.current]))
  end
end

@inline function perform_step!(integrator,cache::StochasticCompositeCache,f::F=integrator.f) where F
  if cache.current == 1
    perform_step!(integrator, @inbounds(cache.caches[1]), f)
  elseif cache.current == 2
    perform_step!(integrator, @inbounds(cache.caches[2]), f)
  else
    perform_step!(integrator, @inbounds(cache.caches[cache.current]), f)
  end
end

@inline choose_algorithm!(integrator,cache::StochasticDiffEqCache) = nothing
@inline function choose_algorithm!(integrator,cache::StochasticCompositeCache)
  new_current = cache.choice_function(integrator)
  if new_current != cache.current
    if new_current == 1
      initialize!(integrator, @inbounds(cache.caches[1]))
    elseif new_current == 2
      initialize!(integrator, @inbounds(cache.caches[2]))
    else
      initialize!(integrator, @inbounds(cache.caches[new_current]))
    end
    if cache.current == 1 && new_current == 2
      reset_alg_dependent_opts!(integrator,integrator.alg.algs[1],integrator.alg.algs[2])
      transfer_cache!(integrator,integrator.cache.caches[1],integrator.cache.caches[2])
    elseif cache.current == 2 && new_current == 1
      reset_alg_dependent_opts!(integrator,integrator.alg.algs[2],integrator.alg.algs[1])
      transfer_cache!(integrator,integrator.cache.caches[2],integrator.cache.caches[1])
    else
      reset_alg_dependent_opts!(integrator,integrator.alg.algs[cache.current],integrator.alg.algs[new_current])
      transfer_cache!(integrator,integrator.cache.caches[cache.current],integrator.cache.caches[new_current])
    end
    cache.current = new_current
  end
end

"""
If no user default, then this will change the default to the defaults
for the second algorithm.
"""
@inline function reset_alg_dependent_opts!(integrator,alg1,alg2)
  integrator.dtchangeable = isdtchangeable(alg2)
  if integrator.opts.adaptive == isadaptive(alg1)
    integrator.opts.adaptive = isadaptive(alg2)
  end
  if integrator.opts.qmin == qmin_default(alg1)
    integrator.opts.qmin = qmin_default(alg2)
  end
  if integrator.opts.qmax == qmax_default(alg1)
    integrator.opts.qmax == qmax_default(alg2)
  end
  if integrator.opts.qsteady_min == qsteady_min_default(alg1)
    integrator.opts.qsteady_min = qsteady_min_default(alg2)
  end
  if integrator.opts.qsteady_max == qsteady_max_default(alg1)
    integrator.opts.qsteady_max = qsteady_max_default(alg2)
  end
  if integrator.opts.controller isa PIController
    if integrator.opts.controller.beta2 == beta2_default(alg1)
      integrator.opts.controller.beta2 = beta2_default(alg2)
    end
    if integrator.opts.controller.beta1 == beta1_default(alg1, integrator.opts.controller.beta2)
      integrator.opts.controller.beta1 = beta1_default(alg2, integrator.opts.controller.beta2)
    end
  end
end

# Write how to transfer the cache variables from one cache to the other
# Example: send the history variables from one multistep method to another

@inline transfer_cache!(integrator,alg1,alg2) = nothing
