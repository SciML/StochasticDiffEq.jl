
struct IIF1MConstantCache{vecuType,rhsType,nl_rhsType} <: StochasticDiffEqConstantCache
  uhold::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

struct IIF1MCache{uType,vecuType,DiffCacheType,rhsType,nl_rhsType,rateType,rateNoiseType,rateNoiseCollectionType,NoiseTmpType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  uhold::vecuType
  dual_cache::DiffCacheType
  tmp::uType
  rhs::rhsType
  nl_rhs::nl_rhsType
  rtmp1::rateType
  rtmp2::rateNoiseType
  rtmp3::rateNoiseCollectionType
  noise_tmp::NoiseTmpType
end

u_cache(c::IIF1MCache)    = (c.uprev2,c.u_old)
du_cache(c::IIF1MCache)   = (c.rtmp1,c.rtmp2,c.rtmp3,c.tmp)
vecu_cache(c::IIF1MCache) = (c.uhold,)
dual_cache(c::IIF1MCache) = (c.dual_cache,)

function alg_cache(alg::IIF1M,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  tmp = zero(rate_prototype)
  rhs = RHS_IIF1M_Scalar(f,tmp,t,t)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  IIF1MConstantCache(uhold,rhs,nl_rhs)
end

function alg_cache(alg::IIF1M,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})

  tmp = similar(u,indices(u)); rtmp1 = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  uhold = vec(u) # this makes uhold the same values as integrator.u
  rhs = RHS_IIF1(f,tmp,t,t,dual_cache,size(u))
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  noise_tmp = tmp

  rtmp2 = zeros(noise_rate_prototype)
  if is_diagonal_noise(prob)
    rtmp3 = rtmp2
  else
    rtmp3 = zeros(rate_prototype)
  end
  IIF1MCache(u,uprev,uhold,dual_cache,tmp,rhs,nl_rhs,rtmp1,rtmp2,rtmp3,noise_tmp)
end

struct IIF2MConstantCache{vecuType,rhsType,nl_rhsType} <: StochasticDiffEqConstantCache
  uhold::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

struct IIF2MCache{uType,vecuType,DiffCacheType,rhsType,nl_rhsType,rateType,rateNoiseType,rateNoiseCollectionType,NoiseTmpType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  uhold::vecuType
  dual_cache::DiffCacheType
  tmp::uType
  rhs::rhsType
  nl_rhs::nl_rhsType
  rtmp1::rateType
  rtmp2::rateNoiseType
  rtmp3::rateNoiseCollectionType
  noise_tmp::NoiseTmpType
end

u_cache(c::IIF2MCache)    = (c.uprev2,c.u_old)
du_cache(c::IIF2MCache)   = (c.rtmp1,c.rtmp2,c.rtmp3,c.tmp)
vecu_cache(c::IIF2MCache) = (c.uhold,)
dual_cache(c::IIF2MCache) = (c.dual_cache,)

function alg_cache(alg::IIF2M,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  tmp = zero(rate_prototype)
  rhs = RHS_IIF2M_Scalar(f,tmp,t,t)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  IIF2MConstantCache(uhold,rhs,nl_rhs)
end

function alg_cache(alg::IIF2M,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})

  tmp = similar(u,indices(u)); rtmp1 = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  uhold = vec(u) # this makes uhold the same values as integrator.u
  rhs = RHS_IIF2(f,tmp,t,t,dual_cache,size(u))
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  noise_tmp = tmp

  rtmp2 = zeros(noise_rate_prototype)
  if is_diagonal_noise(prob)
    rtmp3 = rtmp2
  else
    rtmp3 = zeros(rate_prototype)
  end
  IIF2MCache(u,uprev,uhold,dual_cache,tmp,rhs,nl_rhs,rtmp1,rtmp2,rtmp3,noise_tmp)
end

struct IIF1MilConstantCache{vecuType,rhsType,nl_rhsType} <: StochasticDiffEqConstantCache
  uhold::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end
struct IIF1MilCache{uType,vecuType,DiffCacheType,rhsType,nl_rhsType,rateType,rateNoiseType,rateNoiseCollectionType,NoiseTmpType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  uhold::vecuType
  dual_cache::DiffCacheType
  tmp::uType
  rhs::rhsType
  nl_rhs::nl_rhsType
  rtmp1::rateType
  rtmp2::rateNoiseType
  rtmp3::rateNoiseCollectionType
  noise_tmp::NoiseTmpType
  gtmp::rateNoiseType
  gtmp2::rateNoiseType
end

u_cache(c::IIF1MilCache)    = (c.uprev2,c.u_old)
du_cache(c::IIF1MilCache)   = (c.rtmp1,c.rtmp2,c.rtmp3,c.tmp)
vecu_cache(c::IIF1MilCache) = (c.uhold,)
dual_cache(c::IIF1MilCache) = (c.dual_cache,)

function alg_cache(alg::IIF1Mil,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  C = zero(rate_prototype)
  rhs = RHS_IIF1M_Scalar(f,C,t,t)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  IIF1MilConstantCache(uhold,rhs,nl_rhs)
end

function alg_cache(alg::IIF1Mil,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  tmp = similar(u,indices(u)); rtmp1 = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  uhold = vec(u) # this makes uhold the same values as integrator.u
  rhs = RHS_IIF1(f,tmp,t,t,dual_cache,size(u))
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  noise_tmp = similar(noise_rate_prototype)
  gtmp = similar(noise_rate_prototype); gtmp2 = similar(noise_rate_prototype)
  rtmp2 = zeros(noise_rate_prototype)
  if is_diagonal_noise(prob)
    rtmp3 = rtmp2
  else
    rtmp3 = zeros(rate_prototype)
  end
  IIF1MilCache(u,uprev,uhold,dual_cache,tmp,rhs,nl_rhs,rtmp1,rtmp2,rtmp3,noise_tmp,gtmp,gtmp2)
end
