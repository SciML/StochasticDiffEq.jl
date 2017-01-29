abstract StochasticDiffEqCache <: DECache
abstract StochasticDiffEqConstantCache <: StochasticDiffEqCache
abstract StochasticDiffEqMutableCache <: StochasticDiffEqCache

type StochasticCompositeCache{T,F} <: StochasticDiffEqCache
  caches::T
  choice_function::F
  current::Int
end

function alg_cache{T,algType<:StochasticCompositeAlgorithm}(alg::algType,u,ΔW,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{T}})
  caches = map((x)->alg_cache(x,u,ΔW,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,Val{T}),alg.algs)
  StochasticCompositeCache(caches,alg.choice_function,1)
end

immutable EMConstantCache <: StochasticDiffEqConstantCache end
immutable EMCache{uType} <: StochasticDiffEqMutableCache
  utmp1::uType
  utmp2::uType
end

alg_cache(alg::EM,u,ΔW,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = EMConstantCache()

function alg_cache(alg::EM,u,ΔW,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  utmp1 = similar(u); utmp2 = similar(u)
  EMCache(utmp1,utmp2)
end

immutable RKMilConstantCache <: StochasticDiffEqConstantCache end
immutable RKMilCache{uType} <: StochasticDiffEqMutableCache
  du1::uType
  du2::uType
  K::uType
  utilde::uType
  L::uType
end

alg_cache(alg::RKMil,u,ΔW,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = RKMilConstantCache()

function alg_cache(alg::RKMil,u,ΔW,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  du1 = similar(u); du2 = similar(u)
  K = similar(u); utilde = similar(u); L = similar(u)
  RKMilCache(du1,du2,K,utilde,L)
end

immutable SRA1ConstantCache <: StochasticDiffEqConstantCache end
alg_cache(alg::SRA1,u,ΔW,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = SRA1ConstantCache()

immutable SRA1Cache{randType,uType} <: StochasticDiffEqMutableCache
  chi2::randType
  tmp1::uType
  E₁::uType
  E₂::uType
  gt::uType
  k₁::uType
  k₂::uType
  gpdt::uType
  EEsttmp::uType
end

function alg_cache(alg::SRA1,u,ΔW,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  chi2 = similar(ΔW)
  tmp1 = zeros(u)
  E₁ = zeros(u); gt = zeros(u); gpdt = zeros(u)
  E₂ = zeros(u); k₁ = zeros(u); k₂ = zeros(u)
  EEsttmp = zeros(u)
  SRA1Cache(chi2,tmp1,E₁,E₂,gt,k₁,k₂,gpdt,EEsttmp)
end

immutable SRAConstantCache{VType1,VType2,MType,uType} <: StochasticDiffEqConstantCache
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
  SRAConstantCache(c₀,c₁,A₀,B₀,α,β₁,β₂,stages,H0)
end

function alg_cache(alg::SRA,u,ΔW,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  SRAConstantCache(alg.tableau,rate_prototype)
end

immutable SRACache{uType,tabType} <: StochasticDiffEqMutableCache
  H0::Vector{uType}
  A0temp::uType
  B0temp::uType
  ftmp::uType
  gtmp::uType
  chi2::uType
  atemp::uType
  btemp::uType
  E₁::uType
  E₁temp::uType
  E₂::uType
  EEsttmp::uType
  tab::tabType
end

function alg_cache(alg::SRA,u,ΔW,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  H0 = Vector{typeof(u)}(0)
  tab = SRAConstantCache(alg.tableau,rate_prototype)
  for i = 1:tab.stages
    push!(H0,zeros(u))
  end
  A0temp = zeros(u); B0temp = zeros(u)
  ftmp = zeros(u); gtmp = zeros(u); chi2 = zeros(u)
  atemp = zeros(u); btemp = zeros(u); E₂ = zeros(u); E₁temp = zeros(u)
  E₁ = zeros(u)
  EEsttmp = zeros(u)
  SRACache(H0,A0temp,B0temp,ftmp,gtmp,chi2,atemp,btemp,E₁,E₁temp,E₂,EEsttmp,tab)
end

immutable SRIConstantCache{VType1,VType2,MType,uType} <: StochasticDiffEqConstantCache
  c₀::VType1
  c₁::VType1
  A₀::MType
  A₁::MType
  B₀::MType
  B₁::MType
  α::VType2
  β₁::VType2
  β₂::VType2
  β₃::VType2
  β₄::VType2
  stages::Int
  H0::Vector{uType}
  H1::Vector{uType}
end

function SRIConstantCache(tableau,rate_prototype)
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄ = tableau
  stages = length(α)
  H0 = Array{typeof(rate_prototype)}(stages)
  H1 = Array{typeof(rate_prototype)}(stages)
  SRIConstantCache(c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄,stages,H0,H1)
end

function alg_cache(alg::SRI,u,ΔW,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  SRIConstantCache(alg.tableau,rate_prototype)
end

immutable SRICache{randType,uType,tabType} <: StochasticDiffEqMutableCache
  H0::Vector{uType}
  H1::Vector{uType}
  A0temp::uType
  A1temp::uType
  B0temp::uType
  B1temp::uType
  A0temp2::uType
  A1temp2::uType
  B0temp2::uType
  B1temp2::uType
  atemp::uType
  btemp::uType
  E₁::uType
  E₂::uType
  E₁temp::uType
  ftemp::uType
  gtemp::uType
  chi1::randType
  chi2::randType
  chi3::randType
  EEsttmp::uType
  tab::tabType
end

function alg_cache(alg::SRI,u,ΔW,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  H0 = Vector{typeof(u)}(0)
  H1 = Vector{typeof(u)}(0)
  tab = SRIConstantCache(alg.tableau,rate_prototype)
  for i = 1:tab.stages
    push!(H0,zeros(u))
    push!(H1,zeros(u))
  end
  #TODO Reduce memory
  A0temp = zeros(u); A1temp = zeros(u)
  B0temp = zeros(u); B1temp = zeros(u)
  A0temp2 = zeros(u); A1temp2 = zeros(u)
  B0temp2 = zeros(u); B1temp2 = zeros(u)
  atemp = zeros(u); btemp = zeros(u)
  E₁ = zeros(u); E₂ = zeros(u); E₁temp = zeros(u)
  ftemp = zeros(u); gtemp = zeros(u)
  chi1 = similar(ΔW); chi2 = similar(ΔW); chi3 = similar(ΔW)
  EEsttmp = zeros(u)
  SRICache(H0,H1,A0temp,A1temp,B0temp,B1temp,A0temp2,A1temp2,B0temp2,B1temp2,
    atemp,btemp,E₁,E₂,E₁temp,ftemp,gtemp,chi1,chi2,chi3,EEsttmp,tab)
end

immutable SRIW1ConstantCache <: StochasticDiffEqConstantCache end
alg_cache(alg::SRIW1,u,ΔW,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = SRIW1ConstantCache()

immutable SRIW1Cache{randType,uType} <: StochasticDiffEqMutableCache
  chi1::randType
  chi2::randType
  chi3::randType
  fH01o4::uType
  g₁o2::uType
  H0::uType
  H11::uType
  H12::uType
  H13::uType
  g₂o3::uType
  Fg₂o3::uType
  g₃o3::uType
  Tg₃o3::uType
  mg₁::uType
  E₁::uType
  E₂::uType
  fH01::uType
  fH02::uType
  g₁::uType
  g₂::uType
  g₃::uType
  g₄::uType
  EEsttmp::uType
end

function alg_cache(alg::SRIW1,u,ΔW,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  chi1 = similar(ΔW)
  chi2 = similar(ΔW)
  chi3 = similar(ΔW)
  fH01o4 = zeros(uprev)
  g₁o2 = zeros(u)
  H0 = zeros(u)
  H11 = zeros(u)
  H12 = zeros(u)
  H13 = zeros(u)
  g₂o3 = zeros(u)
  Fg₂o3 = zeros(u)
  g₃o3 = zeros(u)
  Tg₃o3 = zeros(u)
  mg₁ = zeros(u)
  E₁ = zeros(u)
  E₂ = zeros(u)
  fH01 = zeros(u); fH02 = zeros(u)
  g₁ = zeros(u); g₂ = zeros(u); g₃ = zeros(u); g₄ = zeros(u)
  EEsttmp = zeros(u)
  SRIW1Cache(chi1,chi2,chi3,fH01o4,g₁o2,H0,H11,H12,H13,g₂o3,Fg₂o3,g₃o3,Tg₃o3,mg₁,E₁,E₂,fH01,fH02,g₁,g₂,g₃,g₄,EEsttmp)
end
