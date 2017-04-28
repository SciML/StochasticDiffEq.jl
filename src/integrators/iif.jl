type RHS_IIF1M_Scalar{F,CType,tType} <: Function
  f::F
  t::tType
  dt::tType
  tmp::CType
end

function (p::RHS_IIF1M_Scalar)(u,resid)
  resid[1] = u[1] - p.tmp - p.dt*p.f[2](p.t+p.dt,u[1])[1]
end

@inline function initialize!(integrator,cache::Union{IIF1MConstantCache,IIF1MilConstantCache},f=integrator.f)
  cache.uhold[1] = integrator.uprev
end

@inline function perform_step!(integrator,cache::Union{IIF1MConstantCache,IIF1MilConstantCache},f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  @unpack uhold,rhs,nl_rhs = cache
  A = integrator.f[1](t,u)
  if typeof(cache) <: IIF1MilConstantCache
    tmp = expm(A*dt)*(uprev + integrator.g(t,uprev).*W.dW + integrator.g(t,uprev)*(0.87/2)*(W.dW^2 - dt))
  else
    tmp = expm(A*dt)*(uprev + integrator.g(t,uprev).*W.dW)
  end

  if integrator.iter > 1 && !integrator.u_modified
    uhold[1] = current_extrapolant(t+dt,integrator)
  end # else uhold is previous value.

  rhs.t = t
  rhs.dt = dt
  rhs.tmp = tmp
  nlres = integrator.alg.nlsolve(nl_rhs,uhold)

  u = nlres[1]
  @pack integrator = t,dt,u
end

type RHS_IIF1{F,uType,tType,DiffCacheType,SizeType,uidxType} <: Function
  f::F
  tmp::uType
  t::tType
  dt::tType
  dual_cache::DiffCacheType
  sizeu::SizeType
  uidx::uidxType
end
function (p::RHS_IIF1)(u,resid)
  du = get_du(p.dual_cache, eltype(u))
  p.f[2](p.t+p.dt,reshape(u,p.sizeu),du)
  for i in p.uidx
    resid[i] = u[i] - p.tmp[i] - p.dt*du[i]
  end
end

@inline function perform_step!(integrator,cache::Union{IIF1MCache,IIF1MilCache},f=integrator.f)
  @unpack rtmp1,rtmp2,rtmp3,tmp,noise_tmp = cache
  @unpack uhold,rhs,nl_rhs = cache
  @unpack t,dt,uprev,u,W = integrator

  integrator.g(t,uprev,rtmp2)
  if typeof(cache) <: IIF1MCache
    if is_diagonal_noise(integrator.sol.prob)
      for i in eachindex(u)
        rtmp2[i]*=W.dW[i] # rtmp2 === rtmp3
      end
    else
      A_mul_B!(rtmp3,rtmp2,W.dW)
    end
  else #Milstein correction
    for i in eachindex(u)
      noise_tmp[i] = W.dW[i] + (0.87)*(W.dW[i]^2 - dt)/2
    end
    if is_diagonal_noise(integrator.sol.prob)
      for i in eachindex(u)
        rtmp2[i]*= noise_tmp[i] # rtmp2 === rtmp3
      end
    else
      A_mul_B!(rtmp3,rtmp2,noise_tmp)
    end
  end

  for i in eachindex(u)
    rtmp3[i] += uprev[i]
  end

  A = integrator.f[1](t,uprev,rtmp1)
  M = expm(A*dt)
  A_mul_B!(tmp,M,rtmp3)

  rhs.t = t
  rhs.dt = dt
  rhs.tmp = tmp
  rhs.uidx = eachindex(u)
  rhs.sizeu = size(u)
  nlres = integrator.alg.nlsolve(nl_rhs,uhold)

  copy!(uhold,nlres)


  @pack integrator = t,dt,u
end
