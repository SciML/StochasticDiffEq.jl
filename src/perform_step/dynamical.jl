function verify_f2(f, p, q, pa, t, integrator, ::StochasticDynamicalEqConstantCache)
    res = f(p, q, pa, t)
    res != p && throwex(integrator)
end
function verify_f2(f, res, p, q, pa, t, integrator, ::StochasticDiffEqMutableCache)
    f(res, p, q, pa, t)
    res != p && throwex(integrator)
end
function throwex(integrator)
  algn = typeof(integrator.alg)
  throw(ArgumentError("Algorithm $algn is not applicable if f2(p, q, t) != p"))
end

function initialize!(integrator, cache::Union{BAOABConstantCache,ABOBAConstantCache})
  @unpack t,dt,uprev,u,p,W = integrator
  du1 = integrator.uprev.x[1]
  u1 = integrator.uprev.x[2]

  verify_f2(integrator.f.f2, du1, u1, p, t, integrator, cache)
  cache.k = integrator.f.f1(du1,u1,p,t)
end

function initialize!(integrator, cache::Union{BAOABCache,ABOBACache})
  @unpack t,dt,uprev,u,p,W = integrator
  du1 = integrator.uprev.x[1]
  u1 = integrator.uprev.x[2]

  verify_f2(integrator.f.f2, cache.k, du1, u1, p, t, integrator, cache)
  integrator.f.f1(cache.k,du1,u1,p,t)
end

@muladd function perform_step!(integrator,cache::BAOABConstantCache)
  @unpack t,dt,sqdt,uprev,u,p,W,f = integrator
  @unpack half, c1, c2 = cache
  du1 = uprev.x[1]
  u1 = uprev.x[2]

  # B
  du2 = du1 + half*dt*cache.k

  # A
  u2 = u1 + half*dt*du2

  # O
  noise = integrator.g(u2,p,t+dt*half).*W.dW ./ sqdt
  if typeof(c2) <: AbstractMatrix || typeof(noise) <: Number
    du3 = c1*du2 + c2*noise
  else
    du3 = c1.*du2 + c2.*noise
  end

  # A
  u = u2 + half*dt*du3

  # B
  cache.k = f.f1(du3,u,p,t+dt)
  du = du3 + half*dt*cache.k

  integrator.u = ArrayPartition((du, u))
end

@muladd function perform_step!(integrator,cache::BAOABCache)
  @unpack t,dt,sqdt,uprev,u,p,W,f = integrator
  @unpack utmp, dumid, dutmp, dunoise, k, gtmp, noise, half, c1, c2 = cache
  du1 = uprev.x[1]
  u1 = uprev.x[2]

  # B
  @.. dumid = du1 + half*dt*k

  # A
  @.. utmp = u1 + half*dt*dumid

  # O
  integrator.g(gtmp,utmp,p,t+dt*half)
  @.. noise = gtmp*W.dW / sqdt
  if typeof(c2) <: AbstractMatrix
      mul!(dutmp,c1,dumid)
      mul!(dunoise,c2,noise)
      @.. dutmp+= dunoise
  else
      @.. dutmp = c1*dumid + c2*noise
  end

  # A
  @.. u.x[2] = utmp + half*dt*dutmp

  # B
  f.f1(k,dutmp,u.x[2],p,t+dt)
  @.. u.x[1] = dutmp + half*dt*k
end


@muladd function perform_step!(integrator,cache::ABOBAConstantCache)
  @unpack t,dt,sqdt,uprev,u,p,W,f = integrator
  @unpack half, c₂, σ = cache
  du1 = uprev.x[1]
  u1 = uprev.x[2]

  # A
  u_mid = u1 + half*dt*du1

  # BOB: du_t+1
  cache.k = f.f1(du1,u_mid,p,t+half*dt)
  noise = integrator.g(u_mid,p,t+dt*half).*W.dW / sqdt

  if typeof(σ) <: AbstractMatrix || typeof(noise) <: Number
    du = c₂ * (du1 + half*dt .* cache.k) .+ σ*noise .+ half * dt .*cache.k
  else
    du = c₂ .* (du1 + half*dt .* cache.k) .+ σ.*noise .+ half * dt .*cache.k
  end
  # A: xt+1
  u = u_mid .+ half * dt .*du

  integrator.u = ArrayPartition((du, u))
end


@muladd function perform_step!(integrator,cache::ABOBACache)
  @unpack t,dt,sqdt,uprev,u,p,W,f = integrator
  @unpack utmp, dumid, dutmp, dunoise, k, gtmp, noise, half, c₂, σ = cache
  du1 = uprev.x[1]
  u1 = uprev.x[2]

  # A: xt+1/2
  @.. utmp = u1 + half*dt*du1


  # B
  f.f1(k,du1,utmp,p,t+dt)
  @.. dumid = du1 + half*dt*k

  # O
  integrator.g(gtmp,utmp,p,t+dt*half)
  @.. noise = gtmp*W.dW / sqdt

  if typeof(σ) <: AbstractMatrix
      mul!(dutmp,c₂,dumid)
      mul!(dunoise,σ,noise)
      @.. dutmp+=dunoise
  else
      @.. dutmp = c₂*dumid + σ*noise
  end


  # B
  @.. u.x[1] = dutmp + half*dt*k

  # A: xt+1
  @.. u.x[2] = utmp + half*dt*u.x[1]
end




function initialize!(integrator, cache::OBABOConstantCache)
  @unpack t,dt,uprev,u,p,W = integrator
  du1 = integrator.uprev.x[1]
  u1 = integrator.uprev.x[2]

  verify_f2(integrator.f.f2, du1, u1, p, t, integrator, cache)
  cache.k = integrator.f.f1(du1,u1,p,t)
  cache.gt = integrator.g(u1,p,t)
end

function initialize!(integrator, cache::OBABOCache)
  @unpack t,dt,uprev,u,p,W = integrator
  du1 = integrator.uprev.x[1]
  u1 = integrator.uprev.x[2]

  verify_f2(integrator.f.f2, cache.k, du1, u1, p, t, integrator, cache)
  integrator.f.f1(cache.k,du1,u1,p,t)
  integrator.g(cache.gtmp,u1,p,t)
end


@muladd function perform_step!(integrator,cache::OBABOConstantCache)
  @unpack t,dt,sqdt,uprev,u,p,W,f = integrator
  @unpack half, c₂, σ = cache
  du1 = uprev.x[1]
  u1 = uprev.x[2]

  # O
  noise = cache.gt.*W.dW ./ sqdt
  if typeof(σ) <: AbstractMatrix || typeof(noise) <: Number
    du2 = c₂*du1 + σ*noise
  else
    du2 = c₂.*du1 + σ.*noise
  end

  # B
  dumid = du2 + half*dt*cache.k

  # A
  u = u1 + dt*dumid

  cache.k = f.f1(dumid,u,p,t+dt)
  # B
  du3 = dumid + half*dt*cache.k

  # O
  cache.gt = integrator.g(u,p,t+dt)
  noise = cache.gt.*W.dZ ./ sqdt # That should be a second noise
  if typeof(σ) <: AbstractMatrix || typeof(noise) <: Number
    du = c₂*du3 + σ*noise
  else
    du = c₂.*du3 + σ.*noise
  end

  integrator.u = ArrayPartition((du, u))
end


@muladd function perform_step!(integrator,cache::OBABOCache)
  @unpack t,dt,sqdt,uprev,u,p,W,f = integrator
  @unpack utmp, dumid, dutmp, dunoise, k, gtmp, noise, half, c₂, σ = cache
  du1 = uprev.x[1]
  u1 = uprev.x[2]

  # O
  @.. noise = gtmp*W.dW / sqdt

  if typeof(σ) <: AbstractMatrix
      mul!(dutmp,c₂,du1)
      mul!(dunoise,σ,noise)
      @.. dutmp+=dunoise
  else
      @.. dutmp = c₂*du1 + σ*noise
  end

  # B

  @.. dumid = dutmp + half*dt*k

  # A: xt+1
  @.. u.x[2] = u1 + dt*dumid


  # B
  f.f1(k,dumid,u.x[2],p,t+dt)
  @..  dutmp = dumid + half*dt*k

  # O
  integrator.g(gtmp,u.x[2],p,t+dt)
  @.. noise = gtmp*W.dZ / sqdt  # That should be a second noise

  if typeof(σ) <: AbstractMatrix
      mul!(u.x[1],c₂,dutmp)
      mul!(dunoise,σ,noise)
      @.. u.x[1]+=dunoise
  else
      @.. u.x[1] = c₂*dutmp + σ*noise
  end


end
