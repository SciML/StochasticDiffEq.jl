@muladd function perform_step!(integrator,cache::PCEulerConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p,f,g = integrator
  @unpack theta,eta,bbprime = integrator.alg
  if !is_diagonal_noise(integrator.sol.prob) || typeof(W.dW) <: Number
    noises_n = g(uprev,p,t)*W.dW
  else
    noises_n = g(uprev,p,t).*W.dW
  end

  f_n = f(uprev,p,t)
  ubar = uprev + dt*f_n + noises_n

  fbar_n = f_n - eta*bbprime(uprev, p, t)

  tnp1 = t + dt
  fbar_np1 = f(ubar, p, tnp1) - eta*bbprime(ubar, p, tnp1)

  if !is_diagonal_noise(integrator.sol.prob) || typeof(W.dW) <: Number
    noises_np1 = g(ubar,p,tnp1)*W.dW
  else
    noises_np1 = g(ubar,p,tnp1).*W.dW
  end

  u .= uprev + (theta*fbar_np1 + (1-theta)*fbar_n)*dt + eta*noises_np1 + (1-eta)*noises_n
end
