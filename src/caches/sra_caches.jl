struct SRA1ConstantCache <: StochasticDiffEqConstantCache end
alg_cache(alg::SRA1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = SRA1ConstantCache()

struct SRA1Cache{randType,rateType,uType,rateNoiseType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  chi2::randType
  tmp1::uType
  E₁::rateType
  E₂::rateType
  gt::rateNoiseType
  k₁::rateType
  k₂::rateType
  gpdt::rateNoiseType
  tmp::uType
end

u_cache(c::SRA1Cache) = ()
du_cache(c::SRA1Cache) = (c.chi2,c.E₁,c.E₂,c.gt,c.k₁,c.k₂,c.gpdt)
user_cache(c::SRA1Cache) = (c.u,c.uprev,c.tmp,c.tmp1)

function alg_cache(alg::SRA1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    chi2 = copy(ΔW)
  else
    chi2 = zero(ΔW)
  end
  tmp1 = zero(u)
  E₁ = zero(rate_prototype); gt = zero(noise_rate_prototype); gpdt = zero(noise_rate_prototype)
  E₂ = zero(rate_prototype); k₁ = zero(rate_prototype); k₂ = zero(rate_prototype)
  tmp = zero(u)
  SRA1Cache(u,uprev,chi2,tmp1,E₁,E₂,gt,k₁,k₂,gpdt,tmp)
end

struct SRA2ConstantCache{T} <: StochasticDiffEqConstantCache
  a21::T
  b21::T
  c02::T
  c11::T
  c12::T
  α1::T
  α2::T
  beta12::T
  beta21::T
  beta22::T
end

function SRA2ConstantCache(uBottomEltype)
  a21 = uBottomEltype(3//4)
  b21 = uBottomEltype(3//2)
  c02 = uBottomEltype(3//4)
  c11 = uBottomEltype(1//3)
  c12 = uBottomEltype(1)
  α1 = uBottomEltype(1//3)
  α2 = uBottomEltype(2//3)
  beta12 = uBottomEltype(1)
  beta21 = uBottomEltype(3//2)
  beta22 = uBottomEltype(-3//2)
  SRA2ConstantCache(a21,b21,c02,c11,c12,α1,α2,beta12,beta21,beta22)
end

function alg_cache(alg::SRA2,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,uEltypeNoUnits,
                   uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  SRA2ConstantCache(uBottomEltype)
end

struct SRA2Cache{uType,randType,tabType,NT,T} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  chi2::randType
  tab::tabType
  g1::NT
  g2::NT
  k1::T
  k2::T
  E₁::T
  E₂::T
  tmp::T
end

function alg_cache(alg::SRA2,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,uEltypeNoUnits,
                   uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    chi2 = copy(ΔW)
  else
    chi2 = zero(ΔW)
  end
  tab = SRA2ConstantCache(uBottomEltype)
  g1 = zero(noise_rate_prototype); g2 = zero(noise_rate_prototype)
  k1 = zero(rate_prototype); k2 = zero(rate_prototype)
  E₁ = zero(rate_prototype); E₂ = zero(rate_prototype)
  tmp = k1
  SRA2Cache(u,uprev,chi2,tab,g1,g2,k1,k2,E₁,E₂,tmp)
end

struct ThreeStageSRAConstantCache{T} <: StochasticDiffEqConstantCache
  a21::T
  a31::T
  a32::T
  b21::T
  b31::T
  b32::T
  c02::T
  c03::T
  c11::T
  c12::T
  c13::T
  α1::T
  α2::T
  α3::T
  beta11::T
  beta12::T
  beta13::T
  beta21::T
  beta22::T
  beta23::T
end

function SRA3ConstantCache(uBottomEltype)
  a21 = uBottomEltype(1)
  a31 = uBottomEltype(1//4)
  a32 = uBottomEltype(1//4)
  b21 = uBottomEltype(0)
  b31 = uBottomEltype(1)
  b32 = uBottomEltype(1//2)
  c02 = uBottomEltype(1)
  c03 = uBottomEltype(1//2)
  c11 = uBottomEltype(1)
  c12 = uBottomEltype(0)
  c13 = uBottomEltype(0)
  α1 = uBottomEltype(1//6)
  α2 = uBottomEltype(1//6)
  α3 = uBottomEltype(2//3)
  beta11 = uBottomEltype(1)
  beta12 = uBottomEltype(0)
  beta13 = uBottomEltype(0)
  beta21 = uBottomEltype(-1)
  beta22 = uBottomEltype(1)
  beta23 = uBottomEltype(0)
  ThreeStageSRAConstantCache(a21,a31,a32,b21,b31,b32,c02,c03,c11,c12,c13,
                             α1,α2,α3,beta11,beta12,beta13,beta21,beta22,beta23)
end

function SOSRAConstantCache(uBottomEltype)
  α1 = uBottomEltype(0.2889874966892885)
  α2 = uBottomEltype(0.6859880440839937)
  α3 = uBottomEltype(0.025024459226717772)
  c02 = uBottomEltype(0.6923962376159507)
  c03 = uBottomEltype(1)
  c11 = uBottomEltype(0)
  c12 = uBottomEltype(0.041248171110700504)
  c13 = uBottomEltype(1)
  beta11 = uBottomEltype(-16.792534242221663)
  beta12 = uBottomEltype(17.514995785380226)
  beta13 = uBottomEltype(0.27753845684143835)
  beta21 = uBottomEltype(0.4237535769069274)
  beta22 = uBottomEltype(0.6010381474428539)
  beta23 = uBottomEltype(-1.0247917243497813)
  a21 = uBottomEltype(0.6923962376159507)
  a31 = uBottomEltype(-3.1609142252828395)
  a32 = uBottomEltype(4.1609142252828395)
  b21 = uBottomEltype(1.3371632704399763)
  b31 = uBottomEltype(1.442371048468624)
  b32 = uBottomEltype(1.8632741501139225)
  ThreeStageSRAConstantCache(a21,a31,a32,b21,b31,b32,c02,c03,c11,c12,c13,
                             α1,α2,α3,beta11,beta12,beta13,beta21,beta22,beta23)
end

function SOSRA2ConstantCache(uBottomEltype)
  α1 = uBottomEltype(0.4999999999999998)
  α2 = uBottomEltype(-0.9683897375354181)
  α3 = uBottomEltype(1.4683897375354185)
  c02 = uBottomEltype(1)
  c03 = uBottomEltype(1)
  c11 = uBottomEltype(0)
  c12 = uBottomEltype(1)
  c13 = uBottomEltype(1)
  beta11 = uBottomEltype(0)
  beta12 = uBottomEltype(0.92438032145683)
  beta13 = uBottomEltype(0.07561967854316998)
  beta21 = uBottomEltype(1)
  beta22 = uBottomEltype(-0.8169981105823436)
  beta23 = uBottomEltype(-0.18300188941765633)
  a21 = uBottomEltype(1)
  a31 = uBottomEltype(0.9511849235504364)
  a32 = uBottomEltype(0.04881507644956362)
  b21 = uBottomEltype(0.7686101171003622)
  b31 = uBottomEltype(0.43886792994934987)
  b32 = uBottomEltype(0.7490415909204886)
  ThreeStageSRAConstantCache(a21,a31,a32,b21,b31,b32,c02,c03,c11,c12,c13,
                             α1,α2,α3,beta11,beta12,beta13,beta21,beta22,beta23)
end

function alg_cache(alg::SRA3,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,
                   f,t,::Type{Val{false}})
  SRA3ConstantCache(uBottomEltype)
end

function alg_cache(alg::SOSRA,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,
                   f,t,::Type{Val{false}})
  SOSRAConstantCache(uBottomEltype)
end

function alg_cache(alg::SOSRA2,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,
                   f,t,::Type{Val{false}})
  SOSRA2ConstantCache(uBottomEltype)
end

struct ThreeStageSRACache{uType,randType,tabType,NT,T,GT} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  chi2::randType
  tab::tabType
  g1::NT
  g2::NT
  g3::NT
  k1::T
  k2::T
  k3::T
  E₁::T
  E₂::T
  tmp::T
  gtmp::GT
end

function alg_cache(alg::SRA3,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,uEltypeNoUnits,
                   uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    chi2 = copy(ΔW)
  else
    chi2 = zero(ΔW)
  end
  tab = SRA3ConstantCache(uBottomEltype)
  g1 = zero(noise_rate_prototype); g2 = zero(noise_rate_prototype)
  g3 = zero(noise_rate_prototype)
  k1 = zero(rate_prototype); k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  E₁ = zero(rate_prototype); E₂ = zero(rate_prototype)
  tmp = k1
  if typeof(noise_rate_prototype) == typeof(rate_prototype)
    gtmp = nothing
  else
    gtmp = zero(noise_rate_prototype)
  end
  ThreeStageSRACache(u,uprev,chi2,tab,g1,g2,g3,k1,k2,k3,E₁,E₂,tmp,gtmp)
end

function alg_cache(alg::SOSRA,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,uEltypeNoUnits,
                   uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    chi2 = copy(ΔW)
  else
    chi2 = zero(ΔW)
  end
  tab = SOSRAConstantCache(uBottomEltype)
  g1 = zero(noise_rate_prototype); g2 = zero(noise_rate_prototype)
  g3 = zero(noise_rate_prototype)
  k1 = zero(rate_prototype); k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  E₁ = zero(rate_prototype); E₂ = zero(rate_prototype)
  tmp = k1
  if typeof(noise_rate_prototype) == typeof(rate_prototype)
    gtmp = nothing
  else
    gtmp = zero(noise_rate_prototype)
  end
  ThreeStageSRACache(u,uprev,chi2,tab,g1,g2,g3,k1,k2,k3,E₁,E₂,tmp,gtmp)
end

function alg_cache(alg::SOSRA2,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,uEltypeNoUnits,
                   uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    chi2 = copy(ΔW)
  else
    chi2 = zero(ΔW)
  end
  tab = SOSRA2ConstantCache(uBottomEltype)
  g1 = zero(noise_rate_prototype); g2 = zero(noise_rate_prototype)
  g3 = zero(noise_rate_prototype)
  k1 = zero(rate_prototype); k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  E₁ = zero(rate_prototype); E₂ = zero(rate_prototype)
  tmp = k1
  if typeof(noise_rate_prototype) == typeof(rate_prototype)
    gtmp = nothing
  else
    gtmp = zero(noise_rate_prototype)
  end
  ThreeStageSRACache(u,uprev,chi2,tab,g1,g2,g3,k1,k2,k3,E₁,E₂,tmp,gtmp)
end

struct SRAConstantCache{VType1,VType2,MType,uType} <: StochasticDiffEqConstantCache
  c₀::VType1
  c₁::VType1
  A₀::MType
  B₀::MType
  α::VType2
  β₁::VType2
  β₂::VType2
  stages::Int
  H0::Vector{uType}
end

function SRAConstantCache(tableau,rate_prototype)
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂ = tableau
  stages = length(α)
  H0 = Vector{typeof(rate_prototype)}(undef,stages)
  SRAConstantCache(c₀,c₁,A₀',B₀',α,β₁,β₂,stages,H0)
end

function alg_cache(alg::SRA,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  SRAConstantCache(alg.tableau,rate_prototype)
end

struct SRACache{uType,rateType,tabType,randType,rateNoiseType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  H0::Vector{uType}
  A0temp::rateType
  B0temp::rateType
  ftmp::rateType
  gtmp::rateNoiseType
  chi2::randType
  atemp::rateType
  btemp::rateType
  E₁::rateType
  E₁temp::rateType
  E₂::rateType
  tmp::uType
  tab::tabType
end

u_cache(c::SRACache) = ()
du_cache(c::SRACache) = (c.A0temp,c.B0temp,c.ftmp,c.gtmp,c.chi2,c.chi2,c.atemp,
                          c.btemp,c.E₁,c.E₁temp,c.E₂)
user_cache(c::SRACache) = (c.u,c.uprev,c.tmp,c.H0...)

function alg_cache(alg::SRA,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  H0 = Vector{typeof(u)}()
  tab = SRAConstantCache(alg.tableau,rate_prototype)
  for i = 1:tab.stages
    push!(H0,zero(u))
  end
  A0temp = zero(rate_prototype); B0temp = zero(rate_prototype)
  ftmp = zero(rate_prototype); gtmp = zero(noise_rate_prototype);
  if typeof(ΔW) <: Union{SArray,Number}
    chi2 = copy(ΔW)
  else
    chi2 = zero(ΔW)
  end
  atemp = zero(rate_prototype); btemp = zero(rate_prototype); E₂ = zero(rate_prototype); E₁temp = zero(rate_prototype)
  E₁ = zero(rate_prototype)
  tmp = zero(u)
  SRACache(u,uprev,H0,A0temp,B0temp,ftmp,gtmp,chi2,atemp,btemp,E₁,E₁temp,E₂,tmp,tab)
end
