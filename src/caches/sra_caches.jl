struct SRA1ConstantCache <: StochasticDiffEqConstantCache end
alg_cache(alg::SRA1,prob,u,ΔW,ΔZ,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = SRA1ConstantCache()

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

function alg_cache(alg::SRA1,prob,u,ΔW,ΔZ,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    chi2 = copy(ΔW)
  else
    chi2 = similar(ΔW)
  end
  tmp1 = zeros(u)
  E₁ = zeros(rate_prototype); gt = zeros(noise_rate_prototype); gpdt = zeros(noise_rate_prototype)
  E₂ = zeros(rate_prototype); k₁ = zeros(rate_prototype); k₂ = zeros(rate_prototype)
  tmp = zeros(u)
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

function alg_cache(alg::SRA2,prob,u,ΔW,ΔZ,rate_prototype,
                   noise_rate_prototype,uEltypeNoUnits,
                   uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})

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

function alg_cache(alg::SRA3,prob,u,ΔW,ΔZ,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,
                   f,t,::Type{Val{false}})

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

function alg_cache(alg::SOSRA,prob,u,ΔW,ΔZ,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,
                   f,t,::Type{Val{false}})

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
  H0 = Vector{typeof(rate_prototype)}(stages)
  SRAConstantCache(c₀,c₁,A₀',B₀',α,β₁,β₂,stages,H0)
end

function alg_cache(alg::SRA,prob,u,ΔW,ΔZ,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
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

function alg_cache(alg::SRA,prob,u,ΔW,ΔZ,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  H0 = Vector{typeof(u)}(0)
  tab = SRAConstantCache(alg.tableau,rate_prototype)
  for i = 1:tab.stages
    push!(H0,zeros(u))
  end
  A0temp = zeros(rate_prototype); B0temp = zeros(rate_prototype)
  ftmp = zeros(rate_prototype); gtmp = zeros(noise_rate_prototype);
  if typeof(ΔW) <: Union{SArray,Number}
    chi2 = copy(ΔW)
  else
    chi2 = similar(ΔW)
  end
  atemp = zeros(rate_prototype); btemp = zeros(rate_prototype); E₂ = zeros(rate_prototype); E₁temp = zeros(rate_prototype)
  E₁ = zeros(rate_prototype)
  tmp = zeros(u)
  SRACache(u,uprev,H0,A0temp,B0temp,ftmp,gtmp,chi2,atemp,btemp,E₁,E₁temp,E₂,tmp,tab)
end
