struct SRIConstantCache{VType1,VType2,MType,uType} <: StochasticDiffEqConstantCache
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
  error_terms::Int
end

function SRIConstantCache(tableau,rate_prototype,error_terms)
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄ = tableau
  stages = length(α)
  H0 = Array{typeof(rate_prototype)}(undef,stages)
  H1 = Array{typeof(rate_prototype)}(undef,stages)
  SRIConstantCache(c₀,c₁,A₀',A₁',B₀',B₁',α,β₁,β₂,β₃,β₄,stages,H0,H1,error_terms)
end

function alg_cache(alg::SRI,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  SRIConstantCache(alg.tableau,rate_prototype,alg.error_terms)
end

struct SRICache{randType,uType,rateType,tabType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  H0::Vector{uType}
  H1::Vector{uType}
  A0temp::rateType
  A1temp::rateType
  B0temp::rateType
  B1temp::rateType
  A0temp2::rateType
  A1temp2::rateType
  B0temp2::rateType
  B1temp2::rateType
  atemp::rateType
  btemp::rateType
  E₁::rateType
  E₂::rateType
  E₁temp::rateType
  ftemp::rateType
  gtemp::rateType
  chi1::randType
  chi2::randType
  chi3::randType
  tmp::uType
  tab::tabType
end

u_cache(c::SRICache) = ()
du_cache(c::SRICache) = (c.A0temp,c.A1temp,c.B0temp,c.B1temp,c.A0temp2,c.A1temp2,
                          c.B0temp2,c.B1temp2,c.atemp,c.btemp,c.E₁,c.E₂,c.E₁temp,
                          c.ftemp,c.gtemp,c.chi1,c.chi2,c.chi3)
user_cache(c::SRICache) = (c.u,c.uprev,c.tmp,c.H0...,c.H1...)

function alg_cache(alg::SRI,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  H0 = Vector{typeof(u)}()
  H1 = Vector{typeof(u)}()
  tab = SRIConstantCache(alg.tableau,rate_prototype,alg.error_terms)
  for i = 1:tab.stages
    push!(H0,zero(u))
    push!(H1,zero(u))
  end
  #TODO Reduce memory
  A0temp = zero(rate_prototype); A1temp = zero(rate_prototype)
  B0temp = zero(rate_prototype); B1temp = zero(rate_prototype)
  A0temp2 = zero(rate_prototype); A1temp2 = zero(rate_prototype)
  B0temp2 = zero(rate_prototype); B1temp2 = zero(rate_prototype)
  atemp = zero(rate_prototype); btemp = zero(rate_prototype)
  E₁ = zero(rate_prototype); E₂ = zero(rate_prototype); E₁temp = zero(rate_prototype)
  ftemp = zero(rate_prototype); gtemp = zero(rate_prototype)

  if typeof(ΔW) <: Union{SArray,Number}
    chi1 = copy(ΔW)
    chi2 = copy(ΔW)
    chi3 = copy(ΔW)
  else
    chi1 = zero(ΔW)
    chi2 = zero(ΔW)
    chi3 = zero(ΔW)
  end
  tmp = zero(u)
  SRICache(u,uprev,H0,H1,A0temp,A1temp,B0temp,B1temp,A0temp2,A1temp2,B0temp2,B1temp2,
    atemp,btemp,E₁,E₂,E₁temp,ftemp,gtemp,chi1,chi2,chi3,tmp,tab)
end

struct SRIW1ConstantCache <: StochasticDiffEqConstantCache end
alg_cache(alg::SRIW1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = SRIW1ConstantCache()

struct SRIW1Cache{randType,uType,rateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  chi1::randType
  chi2::randType
  chi3::randType
  fH01o4::rateType
  g₁o2::rateType
  H0::uType
  H11::uType
  H12::uType
  H13::uType
  g₂o3::rateType
  Fg₂o3::rateType
  g₃o3::rateType
  Tg₃o3::rateType
  mg₁::rateType
  E₁::rateType
  E₂::rateType
  fH01::rateType
  fH02::rateType
  g₁::rateType
  g₂::rateType
  g₃::rateType
  g₄::rateType
  tmp::uType
end

u_cache(c::SRIW1Cache) = ()
du_cache(c::SRIW1Cache) = (c.chi1,c.chi2,c.chi3,c.fH01o4,c.g₁o2,c.g₂o3,c.Fg₂o3,c.g₃o3,c.Tg₃o3,c.mg₁,
                          c.E₁,c.E₂,c.fH01,c.fH02,c.g₁,c.g₂,c.g₃,c.g₄)
user_cache(c::SRIW1Cache) = (c.u,c.uprev,c.tmp,c.H0,c.H11,c.H12,c.H13)

function alg_cache(alg::SRIW1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    chi1 = copy(ΔW)
    chi2 = copy(ΔW)
    chi3 = copy(ΔW)
  else
    chi1 = zero(ΔW)
    chi2 = zero(ΔW)
    chi3 = zero(ΔW)
  end
  fH01o4 = zero(rate_prototype)
  g₁o2 = zero(rate_prototype)
  H0 = zero(u)
  H11 = zero(u)
  H12 = zero(u)
  H13 = zero(u)
  g₂o3 = zero(rate_prototype)
  Fg₂o3 = zero(rate_prototype)
  g₃o3 = zero(rate_prototype)
  Tg₃o3 = zero(rate_prototype)
  mg₁ = zero(rate_prototype)
  E₁ = zero(rate_prototype)
  E₂ = zero(rate_prototype)
  fH01 = zero(rate_prototype); fH02 = zero(rate_prototype)
  g₁ = zero(rate_prototype); g₂ = zero(rate_prototype); g₃ = zero(rate_prototype); g₄ = zero(rate_prototype)
  tmp = zero(u)
  SRIW1Cache(u,uprev,chi1,chi2,chi3,fH01o4,g₁o2,H0,H11,H12,H13,g₂o3,Fg₂o3,g₃o3,Tg₃o3,mg₁,E₁,E₂,fH01,fH02,g₁,g₂,g₃,g₄,tmp)
end

struct FourStageSRIConstantCache{T} <: StochasticDiffEqConstantCache
  a021::T
  a031::T
  a032::T
  a041::T
  a042::T
  a043::T
  a121::T
  a131::T
  a132::T
  a141::T
  a142::T
  a143::T
  b021::T
  b031::T
  b032::T
  b041::T
  b042::T
  b043::T
  b121::T
  b131::T
  b132::T
  b141::T
  b142::T
  b143::T
  α1::T
  α2::T
  α3::T
  α4::T
  c02::T
  c03::T
  c04::T
  c11::T
  c12::T
  c13::T
  c14::T
  beta11::T
  beta12::T
  beta13::T
  beta14::T
  beta21::T
  beta22::T
  beta23::T
  beta24::T
  beta31::T
  beta32::T
  beta33::T
  beta34::T
  beta41::T
  beta42::T
  beta43::T
  beta44::T
end

function SRIW2ConstantCache(uBottomEltype)
  a021 = uBottomEltype(1)
  a031 = uBottomEltype(1//4)
  a032 = uBottomEltype(1//4)
  a041 = uBottomEltype(0)
  a042 = uBottomEltype(0)
  a043 = uBottomEltype(0)
  a121 = uBottomEltype(1//4)
  a131 = uBottomEltype(1)
  a132 = uBottomEltype(0)
  a141 = uBottomEltype(0)
  a142 = uBottomEltype(0)
  a143 = uBottomEltype(1//4)
  b021 = uBottomEltype(0)
  b031 = uBottomEltype(1)
  b032 = uBottomEltype(1//2)
  b041 = uBottomEltype(0)
  b042 = uBottomEltype(0)
  b043 = uBottomEltype(0)
  b121 = uBottomEltype(-1//2)
  b131 = uBottomEltype(1)
  b132 = uBottomEltype(0)
  b141 = uBottomEltype(2)
  b142 = uBottomEltype(-1)
  b143 = uBottomEltype(1//2)
  α1 = uBottomEltype(1//6)
  α2 = uBottomEltype(1//6)
  α3 = uBottomEltype(2//3)
  α4 = uBottomEltype(0)
  c02 = uBottomEltype(1)
  c03 = uBottomEltype(1//2)
  c04 = uBottomEltype(0)
  c11 = uBottomEltype(0)
  c12 = uBottomEltype(1//4)
  c13 = uBottomEltype(1)
  c14 = uBottomEltype(1//4)
  beta11 = uBottomEltype(-1)
  beta12 = uBottomEltype(4//3)
  beta13 = uBottomEltype(2//3)
  beta14 = uBottomEltype(0)
  beta21 = uBottomEltype(1)
  beta22 = uBottomEltype(-4//3)
  beta23 = uBottomEltype(1//3)
  beta24 = uBottomEltype(0)
  beta31 = uBottomEltype(2)
  beta32 = uBottomEltype(-4//3)
  beta33 = uBottomEltype(-2//3)
  beta34 = uBottomEltype(0)
  beta41 = uBottomEltype(-2)
  beta42 = uBottomEltype(5//3)
  beta43 = uBottomEltype(-2//3)
  beta44 = uBottomEltype(1)
  FourStageSRIConstantCache(a021,a031,a032,a041,a042,a043,a121,a131,a132,a141,
  a142,a143,b021,b031,b032,b041,b042,b043,b121,b131,b132,b141,b142,b143,α1,α2,
  α3,α4,c02,c03,c04,c11,c12,c13,c14,beta11,beta12,beta13,beta14,beta21,beta22,
  beta23,beta24,beta31,beta32,beta33,beta34,beta41,beta42,beta43,beta44)
end

function alg_cache(alg::SRIW2,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,uEltypeNoUnits,
                   uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  SRIW2ConstantCache(uBottomEltype)
end

function SOSRIConstantCache(uBottomEltype)
  a021 = uBottomEltype(-0.04199224421316468)
  a031 = uBottomEltype(2.842612915017106)
  a032 = uBottomEltype(-2.0527723684000727)
  a041 = uBottomEltype(4.338237071435815)
  a042 = uBottomEltype(-2.8895936137439793)
  a043 = uBottomEltype(2.3017575594644466)
  a121 = uBottomEltype(0.26204282091330466)
  a131 = uBottomEltype(0.20903646383505375)
  a132 = uBottomEltype(-0.1502377115150361)
  a141 = uBottomEltype(0.05836595312746999)
  a142 = uBottomEltype(0.6149440396332373)
  a143 = uBottomEltype(0.08535117634046772)
  b021 = uBottomEltype(-0.21641093549612528)
  b031 = uBottomEltype(1.5336352863679572)
  b032 = uBottomEltype(0.26066223492647056)
  b041 = uBottomEltype(-1.0536037558179159)
  b042 = uBottomEltype(1.7015284721089472)
  b043 = uBottomEltype(-0.20725685784180017)
  b121 = uBottomEltype(-0.5119011827621657)
  b131 = uBottomEltype(2.67767339866713)
  b132 = uBottomEltype(-4.9395031322250995)
  b141 = uBottomEltype(0.15580956238299215)
  b142 = uBottomEltype(3.2361551006624674)
  b143 = uBottomEltype(-1.4223118283355949)
  α1 = uBottomEltype(1.140099274172029)
  α2 = uBottomEltype(-0.6401334255743456)
  α3 = uBottomEltype(0.4736296532772559)
  α4 = uBottomEltype(0.026404498125060714)
  c02 = uBottomEltype(-0.04199224421316468)
  c03 = uBottomEltype(0.7898405466170333)
  c04 = uBottomEltype(3.7504010171562823)
  c11 = uBottomEltype(0.0)
  c12 = uBottomEltype(0.26204282091330466)
  c13 = uBottomEltype(0.05879875232001766)
  c14 = uBottomEltype(0.758661169101175)
  beta11 = uBottomEltype(-1.8453464565104432)
  beta12 = uBottomEltype(2.688764531100726)
  beta13 = uBottomEltype(-0.2523866501071323)
  beta14 = uBottomEltype(0.40896857551684956)
  beta21 = uBottomEltype(0.4969658141589478)
  beta22 = uBottomEltype(-0.5771202869753592)
  beta23 = uBottomEltype(-0.12919702470322217)
  beta24 = uBottomEltype(0.2093514975196336)
  beta31 = uBottomEltype(2.8453464565104425)
  beta32 = uBottomEltype(-2.688764531100725)
  beta33 = uBottomEltype(0.2523866501071322)
  beta34 = uBottomEltype(-0.40896857551684945)
  beta41 = uBottomEltype(0.11522663875443433)
  beta42 = uBottomEltype(-0.57877086147738)
  beta43 = uBottomEltype(0.2857851028163886)
  beta44 = uBottomEltype(0.17775911990655704)
  FourStageSRIConstantCache(a021,a031,a032,a041,a042,a043,a121,a131,a132,a141,
  a142,a143,b021,b031,b032,b041,b042,b043,b121,b131,b132,b141,b142,b143,α1,α2,
  α3,α4,c02,c03,c04,c11,c12,c13,c14,beta11,beta12,beta13,beta14,beta21,beta22,
  beta23,beta24,beta31,beta32,beta33,beta34,beta41,beta42,beta43,beta44)
end

function alg_cache(alg::SOSRI,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,uEltypeNoUnits,
                   uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  SOSRIConstantCache(uBottomEltype)
end

function SOSRI2ConstantCache(uBottomEltype)
  a021 = uBottomEltype(0.13804532298278663)
  a031 = uBottomEltype(0.5818361298250374)
  a032 = uBottomEltype(0.4181638701749618)
  a041 = uBottomEltype(0.4670018408674211)
  a042 = uBottomEltype(0.8046204792187386)
  a043 = uBottomEltype(-0.27162232008616016)
  a121 = uBottomEltype(0.45605532163856893)
  a131 = uBottomEltype(0.7555807846451692)
  a132 = uBottomEltype(0.24441921535482677)
  a141 = uBottomEltype(0.6981181143266059)
  a142 = uBottomEltype(0.3453277086024727)
  a143 = uBottomEltype(-0.04344582292908241)
  b021 = uBottomEltype(0.08852381537667678)
  b031 = uBottomEltype(1.0317752458971061)
  b032 = uBottomEltype(0.4563552922077882)
  b041 = uBottomEltype(1.73078280444124)
  b042 = uBottomEltype(-0.46089678470929774)
  b043 = uBottomEltype(-0.9637509618944188)
  b121 = uBottomEltype(0.6753186815412179)
  b131 = uBottomEltype(-0.07452812525785148)
  b132 = uBottomEltype(-0.49783736486149366)
  b141 = uBottomEltype(-0.5591906709928903)
  b142 = uBottomEltype(0.022696571806569924)
  b143 = uBottomEltype(-0.8984927888368557)
  α1 = uBottomEltype(-0.15036858140642623)
  α2 = uBottomEltype(0.7545275856696072)
  α3 = uBottomEltype(0.686995463807979)
  α4 = uBottomEltype(-0.2911544680711602)
  c02 = uBottomEltype(0.13804532298278663)
  c03 = uBottomEltype(0.9999999999999992)
  c04 = uBottomEltype(0.9999999999999994)
  c11 = uBottomEltype(0.0)
  c12 = uBottomEltype(0.45605532163856893)
  c13 = uBottomEltype(0.999999999999996)
  c14 = uBottomEltype(0.9999999999999962)
  beta11 = uBottomEltype(-0.45315689727309133)
  beta12 = uBottomEltype(0.8330937231303951)
  beta13 = uBottomEltype(0.3792843195533544)
  beta14 = uBottomEltype(0.24077885458934192)
  beta21 = uBottomEltype(-0.4994383733810986)
  beta22 = uBottomEltype(0.9181786186154077)
  beta23 = uBottomEltype(-0.25613778661003145)
  beta24 = uBottomEltype(-0.16260245862427797)
  beta31 = uBottomEltype(1.4531568972730915)
  beta32 = uBottomEltype(-0.8330937231303933)
  beta33 = uBottomEltype(-0.3792843195533583)
  beta34 = uBottomEltype(-0.24077885458934023)
  beta41 = uBottomEltype(-0.4976090683622265)
  beta42 = uBottomEltype(0.9148155835648892)
  beta43 = uBottomEltype(-1.4102107084476505)
  beta44 = uBottomEltype(0.9930041932449877)
  FourStageSRIConstantCache(a021,a031,a032,a041,a042,a043,a121,a131,a132,a141,
  a142,a143,b021,b031,b032,b041,b042,b043,b121,b131,b132,b141,b142,b143,α1,α2,
  α3,α4,c02,c03,c04,c11,c12,c13,c14,beta11,beta12,beta13,beta14,beta21,beta22,
  beta23,beta24,beta31,beta32,beta33,beta34,beta41,beta42,beta43,beta44)
end

function alg_cache(alg::SOSRI2,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,uEltypeNoUnits,
                   uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  SOSRI2ConstantCache(uBottomEltype)
end

struct FourStageSRICache{uType,randType,tabType,NT,T} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  chi1::randType
  chi2::randType
  chi3::randType
  tab::tabType
  g1::NT
  g2::NT
  g3::NT
  g4::NT
  k1::T
  k2::T
  k3::T
  k4::T
  E₁::T
  E₂::T
  tmp::T
  H02::T
  H03::T
end

function alg_cache(alg::SRIW2,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,uEltypeNoUnits,
                   uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    chi1 = copy(ΔW)
    chi2 = copy(ΔW)
    chi3 = copy(ΔW)
  else
    chi1 = zero(ΔW)
    chi2 = zero(ΔW)
    chi3 = zero(ΔW)
  end
  tab = SRIW2ConstantCache(uBottomEltype)
  g1 = zero(noise_rate_prototype); g2 = zero(noise_rate_prototype)
  g3 = zero(noise_rate_prototype); g4 = zero(noise_rate_prototype)
  k1 = zero(rate_prototype); k2 = zero(rate_prototype)
  k3 = zero(rate_prototype); k4 = zero(rate_prototype)
  E₁ = zero(rate_prototype); E₂ = zero(rate_prototype)
  tmp = zero(rate_prototype)
  FourStageSRICache(u,uprev,chi1,chi2,chi3,tab,g1,g2,g3,g4,k1,k2,k3,k4,E₁,E₂,tmp,tmp,tmp)
end

function alg_cache(alg::SOSRI,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,uEltypeNoUnits,
                   uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    chi1 = copy(ΔW)
    chi2 = copy(ΔW)
    chi3 = copy(ΔW)
  else
    chi1 = zero(ΔW)
    chi2 = zero(ΔW)
    chi3 = zero(ΔW)
  end
  tab = SOSRIConstantCache(uBottomEltype)
  g1 = zero(noise_rate_prototype); g2 = zero(noise_rate_prototype)
  g3 = zero(noise_rate_prototype); g4 = zero(noise_rate_prototype)
  k1 = zero(rate_prototype); k2 = zero(rate_prototype)
  k3 = zero(rate_prototype); k4 = zero(rate_prototype)
  E₁ = zero(rate_prototype); E₂ = zero(rate_prototype)
  tmp = zero(rate_prototype)
  FourStageSRICache(u,uprev,chi1,chi2,chi3,tab,g1,g2,g3,g4,k1,k2,k3,k4,E₁,E₂,tmp,tmp,tmp)
end

function alg_cache(alg::SOSRI2,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,uEltypeNoUnits,
                   uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    chi1 = copy(ΔW)
    chi2 = copy(ΔW)
    chi3 = copy(ΔW)
  else
    chi1 = zero(ΔW)
    chi2 = zero(ΔW)
    chi3 = zero(ΔW)
  end
  tab = SOSRI2ConstantCache(uBottomEltype)
  g1 = zero(noise_rate_prototype); g2 = zero(noise_rate_prototype)
  g3 = zero(noise_rate_prototype); g4 = zero(noise_rate_prototype)
  k1 = zero(rate_prototype); k2 = zero(rate_prototype)
  k3 = zero(rate_prototype); k4 = zero(rate_prototype)
  E₁ = zero(rate_prototype); E₂ = zero(rate_prototype)
  tmp = zero(rate_prototype); H02 = zero(rate_prototype)
  H03 = zero(rate_prototype)
  FourStageSRICache(u,uprev,chi1,chi2,chi3,tab,g1,g2,g3,g4,k1,k2,k3,k4,E₁,E₂,tmp,H02,H03)
end
