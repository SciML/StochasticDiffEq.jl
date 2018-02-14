@muladd function perform_step!(integrator,cache::PCEulerConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p,f,g = integrator
  @unpack theta,eta,bbprime = integrator.alg
  dW = W.dW
  if !is_diagonal_noise(integrator.sol.prob) || typeof(W.dW) <: Number
    gdW_n = g(uprev,p,t)*dW
  else
    gdW_n = g(uprev,p,t).*dW
  end

  f_n = f(uprev,p,t)
  ubar = @muladd uprev + dt*f_n + gdW_n

  fbar_n = @muladd f_n - eta*bbprime(uprev, p, t)

  tnp1 = t + dt
  fbar_np1 = @muladd f(ubar, p, tnp1) - eta*bbprime(ubar, p, tnp1)

  if !is_diagonal_noise(integrator.sol.prob) || typeof(dW) <: Number
    gdW_np1 = g(ubar,p,tnp1)*dW
  else
    gdW_np1 = g(ubar,p,tnp1).*dW
  end

  u .= @muladd uprev + (theta*dt)*fbar_np1 + ((1-theta)*dt)*fbar_n + eta*gdW_np1 + (1-eta)*gdW_n
end

@muladd function perform_step!(integrator,cache::PCEulerCache,f=integrator.f)
  @unpack utmp,ftmp,gtmp,gdWtmp, bbprimetmp = cache
  @unpack t,dt,uprev,u,W,p,f,g = integrator
  @unpack theta,eta,bbprime = integrator.alg
  dW = W.dW

  f(ftmp,uprev,p,t)
  g(gtmp,uprev,p,t)
  bbprime(bbprimetmp,uprev,p,t)

  if is_diagonal_noise(integrator.sol.prob)
    @tight_loop_macros for i in eachindex(u)
      @inbounds gtmp[i]*=dW[i] # gtmp === gdWtmp
    end
  else
    A_mul_B!(gdWtmp,gtmp,dW)
  end

  @tight_loop_macros for i in eachindex(utmp)
    @inbounds utmp[i] = @muladd uprev[i] + ftmp[i]*dt + gdWtmp[i]
  end

  @tight_loop_macros for i in eachindex(utmp)
    @inbounds u[i] = @muladd uprev[i] + (1-theta)*(ftmp[i]- eta*bbprimetmp[i])*dt + (1-eta)*gdWtmp[i]
  end

  tnp1 = t + dt
  f(ftmp,utmp,p,tnp1)
  g(gtmp,utmp,p,tnp1)
  bbprime(bbprimetmp,utmp,p,tnp1)

  if is_diagonal_noise(integrator.sol.prob)
    @tight_loop_macros for i in eachindex(u)
      @inbounds gtmp[i]*=dW[i] # gtmp === gdWtmp
    end
  else
    A_mul_B!(gdWtmp,gtmp,dW)
  end

  @tight_loop_macros for i in eachindex(utmp)
    @inbounds u[i] += @muladd theta*(ftmp[i]- eta*bbprimetmp[i])*dt + eta*gdWtmp[i]
  end
end
