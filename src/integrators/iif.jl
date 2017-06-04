type RHS_IIF1M_Scalar{F,CType,tType} <: Function
  f::F
  t::tType
  dt::tType
  tmp::CType
end

function (p::RHS_IIF1M_Scalar)(u,resid)
  resid[1] .= u[1] .- p.tmp .- p.dt.*p.f[2](p.t+p.dt,u[1])[1]
end

type RHS_IIF2M_Scalar{F,CType,tType} <: Function
  f::F
  t::tType
  dt::tType
  tmp::CType
end

function (p::RHS_IIF2M_Scalar)(u,resid)
  resid[1] .= u[1] .- p.tmp .- 0.5p.dt.*p.f[2](p.t+p.dt,u[1])[1]
end

@inline function initialize!(integrator,cache::Union{IIF1MConstantCache,IIF2MConstantCache,IIF1MilConstantCache},f=integrator.f)
  cache.uhold[1] = integrator.uprev
end

@inline function perform_step!(integrator,cache::Union{IIF1MConstantCache,IIF2MConstantCache,IIF1MilConstantCache},f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  @unpack uhold,rhs,nl_rhs = cache
  A = integrator.f[1](t,u)
  if typeof(cache) <: IIF1MilConstantCache
    error("Milstein correction does not work.")
  elseif typeof(cache) <: IIF1MConstantCache
    tmp = expm(A*dt)*(uprev .+ integrator.g(t,uprev).*W.dW)
  elseif typeof(cache) <: IIF2MConstantCache
    tmp = expm(A*dt)*(uprev .+ 0.5dt.*integrator.f[2](t,uprev) .+ integrator.g(t,uprev).*W.dW)
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
  @. resid = u - p.tmp - p.dt*du
end

type RHS_IIF2{F,uType,tType,DiffCacheType,SizeType,uidxType} <: Function
  f::F
  tmp::uType
  t::tType
  dt::tType
  dual_cache::DiffCacheType
  sizeu::SizeType
  uidx::uidxType
end
function (p::RHS_IIF2)(u,resid)
  du = get_du(p.dual_cache, eltype(u))
  p.f[2](p.t+p.dt,reshape(u,p.sizeu),du)
  @. resid = u - p.tmp - 0.5p.dt*du
end

@inline function perform_step!(integrator,cache::Union{IIF1MCache,IIF2MCache,IIF1MilCache},f=integrator.f)
  @unpack rtmp1,rtmp2,rtmp3,tmp,noise_tmp = cache
  @unpack uhold,rhs,nl_rhs = cache
  @unpack t,dt,uprev,u,W = integrator

  integrator.g(t,uprev,rtmp2)
  if typeof(cache) <: Union{IIF1MCache,IIF2MCache}
    if is_diagonal_noise(integrator.sol.prob)
      rtmp2 .*=W.dW # rtmp2 === rtmp3
    else
      A_mul_B!(rtmp3,rtmp2,W.dW)
    end
  else #Milstein correction
    error("Milstein correction does not work.")
  end

  rtmp3 .+= uprev

  if typeof(cache) <: IIF2MCache
    integrator.f[2](t,uprev,rtmp1)
    @. rtmp3 = @muladd 0.5dt*rtmp1 + rtmp3
  end

  A = integrator.f[1](t,uprev,rtmp1)
  M = expm(A*dt)
  A_mul_B!(tmp,M,rtmp3)

  if integrator.iter > 1 && !integrator.u_modified
    current_extrapolant!(uhold,t+dt,integrator)
  end # else uhold is previous value.

  rhs.t = t
  rhs.dt = dt
  rhs.tmp = tmp
  rhs.uidx = eachindex(u)
  rhs.sizeu = size(u)
  nlres = integrator.alg.nlsolve(nl_rhs,uhold)

  copy!(uhold,nlres)


  @pack integrator = t,dt,u
end
