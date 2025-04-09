# Perform Step
@muladd function perform_step!(integrator, cache::TauLeapingConstantCache)
  @unpack t, dt, uprev, u, p, P, c = integrator
  
  P === nothing && error("TauLeaping requires a JumpProblem with a RegularJump")
  P.dt = dt
  tmp = c(uprev, p, t, P.dW, nothing)
  integrator.u = uprev .+ tmp

  if integrator.alg.adaptive
    oldrate = P.cache.currate
    newrate = P.cache.rate(integrator.u, p, t + dt)
    EEstcache = @. abs(newrate - oldrate) / max(50 * integrator.opts.reltol * oldrate, integrator.rate_constants / integrator.dt)
    integrator.EEst = integrator.opts.internalnorm(EEstcache, t)
    if integrator.EEst <= 1
      P.cache.currate = newrate
    end
  else
    integrator.EEst = 1.0
  end
end

@muladd function perform_step!(integrator, cache::TauLeapingCache)
  @unpack t, dt, uprev, u, p, P, c = integrator
  @unpack tmp, rate, newrate, EEstcache = cache
  
  P === nothing && error("TauLeaping requires a JumpProblem with a RegularJump")
  P.dt = dt
  c(tmp, uprev, p, t, P.dW, nothing)
  @.. u = uprev + tmp

  if integrator.alg.adaptive
    P.cache.rate(newrate, u, p, t + dt)
    P.cache.rate(rate, uprev, p, t)
    @.. EEstcache = abs(newrate - rate) / max(50 * integrator.opts.reltol * rate, integrator.rate_constants / integrator.dt)
    integrator.EEst = integrator.opts.internalnorm(EEstcache, t)
    if integrator.EEst <= 1
      P.cache.currate .= newrate
    end
  else
    integrator.EEst = 1.0
  end
end

@muladd function perform_step!(integrator, cache::CaoTauLeapingConstantCache)
  @unpack t, dt, uprev, u, p, P, c = integrator
  
  P === nothing && error("CaoTauLeaping requires a JumpProblem with a RegularJump")
  P.dt = dt
  tmp = c(uprev, p, t, P.dW, nothing)
  integrator.u = uprev .+ tmp

  integrator.EEst = 1.0
end

@muladd function perform_step!(integrator, cache::CaoTauLeapingCache)
  @unpack t, dt, uprev, u, p, P, c = integrator
  @unpack tmp = cache
  
  P === nothing && error("CaoTauLeaping requires a JumpProblem with a RegularJump")
  P.dt = dt
  c(tmp, uprev, p, t, P.dW, nothing)
  @.. u = uprev + tmp

  integrator.EEst = 1.0
end
