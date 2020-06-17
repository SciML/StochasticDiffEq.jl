################################################################################
# Roessler SRK for second order weak approx

struct DRI1ConstantCache{T,T2} <: StochasticDiffEqConstantCache
  # hard-coded version
  a021::T
  a031::T
  a032::T

  a121::T
  a131::T
  #a132::T

  #a221::T
  #a231::T
  #a232::T

  b021::T
  b031::T
  #b032::T

  b121::T
  b131::T
  #b132::T

  b221::T
  b222::T
  b223::T
  b231::T
  b232::T
  b233::T

  α1::T
  α2::T
  α3::T

  c02::T2
  c03::T2

  #c11::T2
  c12::T2
  c13::T2

  #c21::T2
  #c22::T2
  #c23::T2

  beta11::T
  beta12::T
  beta13::T

  #beta21::T
  beta22::T
  beta23::T

  beta31::T
  beta32::T
  beta33::T

  #beta41::T
  beta42::T
  beta43::T

  #quantile(Normal(),1/6)
  NORMAL_ONESIX_QUANTILE::T
end

function DRI1ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}

  a021 = convert(T, 1//2)
  a031 = convert(T, -1)
  a032 = convert(T, 2)

  a121 = convert(T, 342//491)
  a131 = convert(T, 342//491)

  b021 = convert(T, (6-sqrt(6))/10)
  b031 = convert(T, (3+2*sqrt(6))/5)

  b121 = convert(T, 3*sqrt(38//491))
  b131 = convert(T, -3*sqrt(38/491))

  b221 = convert(T, -214//513*sqrt(1105//991))
  b222 = convert(T, -491//513*sqrt(221//4955))
  b223 = convert(T, -491//513*sqrt(221//4955))
  b231 = convert(T, 214//513*sqrt(1105//991))
  b232 = convert(T, 491//513*sqrt(221//4955))
  b233 = convert(T, 491//513*sqrt(221//4955))

  α1 = convert(T, 1//6)
  α2 = convert(T, 2//3)
  α3 = convert(T, 1//6)

  c02 = convert(T2, 1//2)
  c03 = convert(T2, 1)

  c12 = convert(T2, 342//491)
  c13 = convert(T2, 342//491)

  beta11 = convert(T, 193//684)
  beta12 = convert(T, 491//1368)
  beta13 = convert(T, 491//1368)

  beta22 = convert(T, 1//6*sqrt(491//38))
  beta23 = convert(T, -1//6*sqrt(491//38))

  beta31 = convert(T, -4955//7072)
  beta32 = convert(T, 4955//14144)
  beta33 = convert(T, 4955//14144)

  beta42 = convert(T, -1//8*sqrt(4955//221))
  beta43 = convert(T, 1//8*sqrt(4955//221))

  NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

  DRI1ConstantCache(a021,a031,a032,a121,a131,b021,b031,b121,b131,b221,b222,b223,b231,b232,b233,α1,α2,α3,c02,c03,c12,c13,beta11,beta12,beta13,beta22,beta23,beta31,beta32,beta33,beta42,beta43,NORMAL_ONESIX_QUANTILE)
end


function alg_cache(alg::DRI1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  DRI1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end



function RI1ConstantCache(T::Type, T2::Type)
  a021 = convert(T, 2//3) # convert(T, 2//3)
  a031 = convert(T, -1//3)
  a032 = convert(T, 1)

  a121 = convert(T, 1)
  a131 = convert(T, 1)

  b021 = convert(T, 1)
  b031 = convert(T, 0)

  b121 = convert(T, 1)
  b131 = convert(T, -1)

  b221 = convert(T, 1)
  b222 = convert(T, 0)
  b223 = convert(T, 0)
  b231 = convert(T, -1)
  b232 = convert(T, 0)
  b233 = convert(T, 0)

  α1 = convert(T, 1//4)
  α2 = convert(T, 1//2)
  α3 = convert(T, 1//4)

  c02 = convert(T2, 2//3)
  c03 = convert(T2, 2//3)

  c12 = convert(T2, 1)
  c13 = convert(T2, 1)

  beta11 = convert(T, 1//2)
  beta12 = convert(T, 1//4)
  beta13 = convert(T, 1//4)

  beta22 = convert(T, 1//2)
  beta23 = convert(T, -1//2)

  beta31 = convert(T, -1//2)
  beta32 = convert(T, 1//4)
  beta33 = convert(T, 1//4)

  beta42 = convert(T, 1//2)
  beta43 = convert(T, -1//2)

  NORMAL_ONESIX_QUANTILE = convert(T,-0.9674215661017014)

  DRI1ConstantCache(a021,a031,a032,a121,a131,b021,b031,b121,b131,b221,b222,b223,b231,b232,b233,α1,α2,α3,c02,c03,c12,c13,beta11,beta12,beta13,beta22,beta23,beta31,beta32,beta33,beta42,beta43,NORMAL_ONESIX_QUANTILE)
end


function alg_cache(alg::RI1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  RI1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function RI3ConstantCache(T::Type, T2::Type)
  a021 = convert(T, 1)
  a031 = convert(T, 1//4)
  a032 = convert(T, 1//4)

  a121 = convert(T, 1)
  a131 = convert(T, 1)

  b021 = convert(T, (3-2*sqrt(6))/5)
  b031 = convert(T, (6+sqrt(6))/10)

  b121 = convert(T, 1)
  b131 = convert(T, -1)

  b221 = convert(T, 1)
  b222 = convert(T, 0)
  b223 = convert(T, 0)
  b231 = convert(T, -1)
  b232 = convert(T, 0)
  b233 = convert(T, 0)

  α1 = convert(T, 1//6)
  α2 = convert(T, 1//6)
  α3 = convert(T, 2//3)

  c02 = convert(T2, 1)
  c03 = convert(T2, 1//2)

  c12 = convert(T2, 1)
  c13 = convert(T2, 1)

  beta11 = convert(T, 1//2)
  beta12 = convert(T, 1//4)
  beta13 = convert(T, 1//4)

  beta22 = convert(T, 1//2)
  beta23 = convert(T, -1//2)

  beta31 = convert(T, -1//2)
  beta32 = convert(T, 1//4)
  beta33 = convert(T, 1//4)

  beta42 = convert(T, 1//2)
  beta43 = convert(T, -1//2)

  NORMAL_ONESIX_QUANTILE = convert(T,-0.9674215661017014)

  DRI1ConstantCache(a021,a031,a032,a121,a131,b021,b031,b121,b131,b221,b222,b223,b231,b232,b233,α1,α2,α3,c02,c03,c12,c13,beta11,beta12,beta13,beta22,beta23,beta31,beta32,beta33,beta42,beta43,NORMAL_ONESIX_QUANTILE)
end


function alg_cache(alg::RI3,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  RI3ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function RI5ConstantCache(T::Type, T2::Type)
  a021 = convert(T, 1)
  a031 = convert(T, 25//144)
  a032 = convert(T, 35//144)

  a121 = convert(T, 1//4)
  a131 = convert(T, 1//4)

  b021 = convert(T, 1//3)
  b031 = convert(T, -5//6)

  b121 = convert(T, 1//2)
  b131 = convert(T, -1//2)

  b221 = convert(T, 1)
  b222 = convert(T, 0)
  b223 = convert(T, 0)
  b231 = convert(T, -1)
  b232 = convert(T, 0)
  b233 = convert(T, 0)

  α1 = convert(T, 1//10)
  α2 = convert(T, 3//14)
  α3 = convert(T, 24//35)

  c02 = convert(T2, 1)
  c03 = convert(T2, 5//12)

  c12 = convert(T2, 1//4)
  c13 = convert(T2, 1//4)

  beta11 = convert(T, 1)
  beta12 = convert(T, -1)
  beta13 = convert(T, -1)

  beta22 = convert(T, 1)
  beta23 = convert(T, -1)

  beta31 = convert(T, 1//2)
  beta32 = convert(T, -1//4)
  beta33 = convert(T, -1//4)

  beta42 = convert(T, 1//2)
  beta43 = convert(T, -1//2)

  NORMAL_ONESIX_QUANTILE = convert(T,-0.9674215661017014)

  DRI1ConstantCache(a021,a031,a032,a121,a131,b021,b031,b121,b131,b221,b222,b223,b231,b232,b233,α1,α2,α3,c02,c03,c12,c13,beta11,beta12,beta13,beta22,beta23,beta31,beta32,beta33,beta42,beta43,NORMAL_ONESIX_QUANTILE)
end


function alg_cache(alg::RI5,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  RI5ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function RI6ConstantCache(T::Type, T2::Type)
  a021 = convert(T, 1)
  a031 = convert(T, 0)
  a032 = convert(T, 0)

  a121 = convert(T, 1)
  a131 = convert(T, 1)

  b021 = convert(T, 1)
  b031 = convert(T, 0)

  b121 = convert(T, 1)
  b131 = convert(T, -1)

  b221 = convert(T, 1)
  b222 = convert(T, 0)
  b223 = convert(T, 0)
  b231 = convert(T, -1)
  b232 = convert(T, 0)
  b233 = convert(T, 0)

  α1 = convert(T, 1//2)
  α2 = convert(T, 1//2)
  α3 = convert(T, 0)

  c02 = convert(T2, 1)
  c03 = convert(T2, 0)

  c12 = convert(T2, 1)
  c13 = convert(T2, 1)

  beta11 = convert(T, 1//2)
  beta12 = convert(T, 1//4)
  beta13 = convert(T, 1//4)

  beta22 = convert(T, 1//2)
  beta23 = convert(T, -1//2)

  beta31 = convert(T, -1//2)
  beta32 = convert(T, 1//4)
  beta33 = convert(T, 1//4)

  beta42 = convert(T, 1//2)
  beta43 = convert(T, -1//2)

  NORMAL_ONESIX_QUANTILE = convert(T,-0.9674215661017014)

  DRI1ConstantCache(a021,a031,a032,a121,a131,b021,b031,b121,b131,b221,b222,b223,b231,b232,b233,α1,α2,α3,c02,c03,c12,c13,beta11,beta12,beta13,beta22,beta23,beta31,beta32,beta33,beta42,beta43,NORMAL_ONESIX_QUANTILE)
end


function alg_cache(alg::RI6,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  RI6ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function RDI2WMConstantCache(T::Type, T2::Type)
  a021 = convert(T, 1)
  a031 = convert(T, 0)
  a032 = convert(T, 0)

  a121 = convert(T, 2//3)
  a131 = convert(T, 2//3)

  b021 = convert(T, 1)
  b031 = convert(T, 0)

  b121 = convert(T, sqrt(2//3))
  b131 = convert(T, -sqrt(2//3))

  b221 = convert(T, sqrt(2))
  b222 = convert(T, 0)
  b223 = convert(T, 0)
  b231 = convert(T, -sqrt(2))
  b232 = convert(T, 0)
  b233 = convert(T, 0)

  α1 = convert(T, 1//2)
  α2 = convert(T, 1//2)
  α3 = convert(T, 0)

  c02 = convert(T2, 1)
  c03 = convert(T2, 0)

  c12 = convert(T2, 2//3)
  c13 = convert(T2, 2//3)

  beta11 = convert(T, 1//4)
  beta12 = convert(T, 3//8)
  beta13 = convert(T, 3//8)

  beta22 = convert(T, sqrt(6)/4)
  beta23 = convert(T, -sqrt(6)/4)

  beta31 = convert(T, -1//4)
  beta32 = convert(T, 1//8)
  beta33 = convert(T, 1//8)

  beta42 = convert(T, sqrt(2)/4)
  beta43 = convert(T, -sqrt(2)/4)

  NORMAL_ONESIX_QUANTILE = convert(T,-0.9674215661017014)

  DRI1ConstantCache(a021,a031,a032,a121,a131,b021,b031,b121,b131,b221,b222,b223,b231,b232,b233,α1,α2,α3,c02,c03,c12,c13,beta11,beta12,beta13,beta22,beta23,beta31,beta32,beta33,beta42,beta43,NORMAL_ONESIX_QUANTILE)
end


function alg_cache(alg::RDI2WM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  RDI2WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function RDI3WMConstantCache(T::Type, T2::Type)
  a021 = convert(T, 1//2)
  a031 = convert(T, 0)
  a032 = convert(T, 3//4)

  a121 = convert(T, 2//3)
  a131 = convert(T, 2//3)

  b021 = convert(T, (9-2*sqrt(15))/14)
  b031 = convert(T, (18+3*sqrt(15))/28)

  b121 = convert(T, sqrt(2//3))
  b131 = convert(T, -sqrt(2//3))

  b221 = convert(T, sqrt(2))
  b222 = convert(T, 0)
  b223 = convert(T, 0)
  b231 = convert(T, -sqrt(2))
  b232 = convert(T, 0)
  b233 = convert(T, 0)

  α1 = convert(T, 2//9)
  α2 = convert(T, 1//3)
  α3 = convert(T, 4//9)

  c02 = convert(T2, 1//2)
  c03 = convert(T2, 3//4)

  c12 = convert(T2, 2//3)
  c13 = convert(T2, 2//3)

  beta11 = convert(T, 1//4)
  beta12 = convert(T, 3//8)
  beta13 = convert(T, 3//8)

  beta22 = convert(T, sqrt(6)/4)
  beta23 = convert(T, -sqrt(6)/4)

  beta31 = convert(T, -1//4)
  beta32 = convert(T, 1//8)
  beta33 = convert(T, 1//8)

  beta42 = convert(T, sqrt(2)/4)
  beta43 = convert(T, -sqrt(2)/4)

  NORMAL_ONESIX_QUANTILE = convert(T,-0.9674215661017014)

  DRI1ConstantCache(a021,a031,a032,a121,a131,b021,b031,b121,b131,b221,b222,b223,b231,b232,b233,α1,α2,α3,c02,c03,c12,c13,beta11,beta12,beta13,beta22,beta23,beta31,beta32,beta33,beta42,beta43,NORMAL_ONESIX_QUANTILE)
end


function alg_cache(alg::RDI3WM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  RDI3WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function RDI4WMConstantCache(T::Type, T2::Type)
  a021 = convert(T, 1//2)
  a031 = convert(T, -1)
  a032 = convert(T, 2)

  a121 = convert(T, 2//3)
  a131 = convert(T, 2//3)

  b021 = convert(T, (6-sqrt(6))/10)
  b031 = convert(T, (3+2*sqrt(6))/5)

  b121 = convert(T, sqrt(2//3))
  b131 = convert(T, -sqrt(2//3))

  b221 = convert(T, sqrt(2))
  b222 = convert(T, 0)
  b223 = convert(T, 0)
  b231 = convert(T, -sqrt(2))
  b232 = convert(T, 0)
  b233 = convert(T, 0)

  α1 = convert(T, 1//6)
  α2 = convert(T, 2//3)
  α3 = convert(T, 1//6)

  c02 = convert(T2, 1//2)
  c03 = convert(T2, 1)

  c12 = convert(T2, 2//3)
  c13 = convert(T2, 2//3)

  beta11 = convert(T, 1//4)
  beta12 = convert(T, 3//8)
  beta13 = convert(T, 3//8)

  beta22 = convert(T, sqrt(6)/4)
  beta23 = convert(T, -sqrt(6)/4)

  beta31 = convert(T, -1//4)
  beta32 = convert(T, 1//8)
  beta33 = convert(T, 1//8)

  beta42 = convert(T, sqrt(2)/4)
  beta43 = convert(T, -sqrt(2)/4)

  NORMAL_ONESIX_QUANTILE = convert(T,-0.9674215661017014)

  DRI1ConstantCache(a021,a031,a032,a121,a131,b021,b031,b121,b131,b221,b222,b223,b231,b232,b233,α1,α2,α3,c02,c03,c12,c13,beta11,beta12,beta13,beta22,beta23,beta31,beta32,beta33,beta42,beta43,NORMAL_ONESIX_QUANTILE)
end


function alg_cache(alg::RDI4WM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  RDI4WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct DRI1Cache{uType,randType,MType1,tabType,rateNoiseType,rateType,possibleRateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  uhat::uType

  _dW::randType
  _dZ::randType
  chi1::randType
  Ihat2::MType1

  tab::tabType

  g1::rateNoiseType
  g2::Vector{rateNoiseType}
  g3::Vector{rateNoiseType}

  k1::rateType
  k2::rateType
  k3::rateType

  H02::uType
  H03::uType
  H12::Vector{uType}
  H13::Vector{uType}
  H22::Vector{uType}
  H23::Vector{uType}

  tmp1::possibleRateType
  tmpg::rateNoiseType

  tmp::uType
  resids::uType

end


function alg_cache(alg::DRI1,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,
                   uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
    _dZ = copy(ΔW)
    chi1 = copy(ΔW)
  else
    _dW = zero(ΔW)
    _dZ = zero(ΔW)
    chi1 = zero(ΔW)
  end
  m = length(ΔW)
  Ihat2 = zeros(eltype(ΔW), m, m)
  tab = DRI1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  g1 = zero(noise_rate_prototype)
  g2 = [zero(noise_rate_prototype) for k=1:m]
  g3 = [zero(noise_rate_prototype) for k=1:m]
  k1 = zero(rate_prototype); k2 = zero(rate_prototype); k3 = zero(rate_prototype)

  H02 = zero(u)
  H03 = zero(u)
  H12 = Vector{typeof(u)}()
  H13 = Vector{typeof(u)}()
  H22 = Vector{typeof(u)}()
  H23 = Vector{typeof(u)}()

  for k=1:m
    push!(H12,zero(u))
    push!(H13,zero(u))
    push!(H22,zero(u))
    push!(H23,zero(u))
  end

  tmp1 = zero(rate_prototype)
  tmpg = zero(noise_rate_prototype)

  uhat = copy(uprev)
  tmp = zero(u)
  resids = zero(u)

  DRI1Cache(u,uprev,uhat,_dW,_dZ,chi1,Ihat2,tab,g1,g2,g3,k1,k2,k3,H02,H03,H12,H13,H22,H23,tmp1,tmpg,tmp,resids)
end


function alg_cache(alg::RI1,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,
                   uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
    _dZ = copy(ΔW)
    chi1 = copy(ΔW)
  else
    _dW = zero(ΔW)
    _dZ = zero(ΔW)
    chi1 = zero(ΔW)
  end
  m = length(ΔW)
  Ihat2 = zeros(eltype(ΔW), m, m)
  tab = RI1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  g1 = zero(noise_rate_prototype)
  g2 = [zero(noise_rate_prototype) for k=1:m]
  g3 = [zero(noise_rate_prototype) for k=1:m]
  k1 = zero(rate_prototype); k2 = zero(rate_prototype); k3 = zero(rate_prototype)

  H02 = zero(u)
  H03 = zero(u)
  H12 = Vector{typeof(u)}()
  H13 = Vector{typeof(u)}()
  H22 = Vector{typeof(u)}()
  H23 = Vector{typeof(u)}()

  for k=1:m
    push!(H12,zero(u))
    push!(H13,zero(u))
    push!(H22,zero(u))
    push!(H23,zero(u))
  end

  tmp1 = zero(rate_prototype)
  tmpg = zero(noise_rate_prototype)

  uhat = copy(uprev)
  tmp = zero(u)
  resids = zero(u)

  DRI1Cache(u,uprev,uhat,_dW,_dZ,chi1,Ihat2,tab,g1,g2,g3,k1,k2,k3,H02,H03,H12,H13,H22,H23,tmp1,tmpg,tmp,resids)
end



function alg_cache(alg::RI3,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,
                   uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
    _dZ = copy(ΔW)
    chi1 = copy(ΔW)
  else
    _dW = zero(ΔW)
    _dZ = zero(ΔW)
    chi1 = zero(ΔW)
  end
  m = length(ΔW)
  Ihat2 = zeros(eltype(ΔW), m, m)
  tab = RI3ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  g1 = zero(noise_rate_prototype)
  g2 = [zero(noise_rate_prototype) for k=1:m]
  g3 = [zero(noise_rate_prototype) for k=1:m]
  k1 = zero(rate_prototype); k2 = zero(rate_prototype); k3 = zero(rate_prototype)

  H02 = zero(u)
  H03 = zero(u)
  H12 = Vector{typeof(u)}()
  H13 = Vector{typeof(u)}()
  H22 = Vector{typeof(u)}()
  H23 = Vector{typeof(u)}()

  for k=1:m
    push!(H12,zero(u))
    push!(H13,zero(u))
    push!(H22,zero(u))
    push!(H23,zero(u))
  end

  tmp1 = zero(rate_prototype)
  tmpg = zero(noise_rate_prototype)

  uhat = copy(uprev)
  tmp = zero(u)
  resids = zero(u)

  DRI1Cache(u,uprev,uhat,_dW,_dZ,chi1,Ihat2,tab,g1,g2,g3,k1,k2,k3,H02,H03,H12,H13,H22,H23,tmp1,tmpg,tmp,resids)
end



function alg_cache(alg::RI5,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,
                   uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
    _dZ = copy(ΔW)
    chi1 = copy(ΔW)
  else
    _dW = zero(ΔW)
    _dZ = zero(ΔW)
    chi1 = zero(ΔW)
  end
  m = length(ΔW)
  Ihat2 = zeros(eltype(ΔW), m, m)
  tab = RI5ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  g1 = zero(noise_rate_prototype)
  g2 = [zero(noise_rate_prototype) for k=1:m]
  g3 = [zero(noise_rate_prototype) for k=1:m]
  k1 = zero(rate_prototype); k2 = zero(rate_prototype); k3 = zero(rate_prototype)

  H02 = zero(u)
  H03 = zero(u)
  H12 = Vector{typeof(u)}()
  H13 = Vector{typeof(u)}()
  H22 = Vector{typeof(u)}()
  H23 = Vector{typeof(u)}()

  for k=1:m
    push!(H12,zero(u))
    push!(H13,zero(u))
    push!(H22,zero(u))
    push!(H23,zero(u))
  end

  tmp1 = zero(rate_prototype)
  tmpg = zero(noise_rate_prototype)

  uhat = copy(uprev)
  tmp = zero(u)
  resids = zero(u)

  DRI1Cache(u,uprev,uhat,_dW,_dZ,chi1,Ihat2,tab,g1,g2,g3,k1,k2,k3,H02,H03,H12,H13,H22,H23,tmp1,tmpg,tmp,resids)
end



function alg_cache(alg::RI6,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,
                   uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
    _dZ = copy(ΔW)
    chi1 = copy(ΔW)
  else
    _dW = zero(ΔW)
    _dZ = zero(ΔW)
    chi1 = zero(ΔW)
  end
  m = length(ΔW)
  Ihat2 = zeros(eltype(ΔW), m, m)
  tab = RI6ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  g1 = zero(noise_rate_prototype)
  g2 = [zero(noise_rate_prototype) for k=1:m]
  g3 = [zero(noise_rate_prototype) for k=1:m]
  k1 = zero(rate_prototype); k2 = zero(rate_prototype); k3 = zero(rate_prototype)

  H02 = zero(u)
  H03 = zero(u)
  H12 = Vector{typeof(u)}()
  H13 = Vector{typeof(u)}()
  H22 = Vector{typeof(u)}()
  H23 = Vector{typeof(u)}()

  for k=1:m
    push!(H12,zero(u))
    push!(H13,zero(u))
    push!(H22,zero(u))
    push!(H23,zero(u))
  end

  tmp1 = zero(rate_prototype)
  tmpg = zero(noise_rate_prototype)

  uhat = copy(uprev)
  tmp = zero(u)
  resids = zero(u)

  DRI1Cache(u,uprev,uhat,_dW,_dZ,chi1,Ihat2,tab,g1,g2,g3,k1,k2,k3,H02,H03,H12,H13,H22,H23,tmp1,tmpg,tmp,resids)
end


function alg_cache(alg::RDI2WM,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,
                   uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
    _dZ = copy(ΔW)
    chi1 = copy(ΔW)
  else
    _dW = zero(ΔW)
    _dZ = zero(ΔW)
    chi1 = zero(ΔW)
  end
  m = length(ΔW)
  Ihat2 = zeros(eltype(ΔW), m, m)
  tab = RDI2WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  g1 = zero(noise_rate_prototype)
  g2 = [zero(noise_rate_prototype) for k=1:m]
  g3 = [zero(noise_rate_prototype) for k=1:m]
  k1 = zero(rate_prototype); k2 = zero(rate_prototype); k3 = zero(rate_prototype)

  H02 = zero(u)
  H03 = zero(u)
  H12 = Vector{typeof(u)}()
  H13 = Vector{typeof(u)}()
  H22 = Vector{typeof(u)}()
  H23 = Vector{typeof(u)}()

  for k=1:m
    push!(H12,zero(u))
    push!(H13,zero(u))
    push!(H22,zero(u))
    push!(H23,zero(u))
  end

  tmp1 = zero(rate_prototype)
  tmpg = zero(noise_rate_prototype)

  uhat = copy(uprev)
  tmp = zero(u)
  resids = zero(u)

  DRI1Cache(u,uprev,uhat,_dW,_dZ,chi1,Ihat2,tab,g1,g2,g3,k1,k2,k3,H02,H03,H12,H13,H22,H23,tmp1,tmpg,tmp,resids)
end




function alg_cache(alg::RDI3WM,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,
                   uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
    _dZ = copy(ΔW)
    chi1 = copy(ΔW)
  else
    _dW = zero(ΔW)
    _dZ = zero(ΔW)
    chi1 = zero(ΔW)
  end
  m = length(ΔW)
  Ihat2 = zeros(eltype(ΔW), m, m)
  tab = RDI3WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  g1 = zero(noise_rate_prototype)
  g2 = [zero(noise_rate_prototype) for k=1:m]
  g3 = [zero(noise_rate_prototype) for k=1:m]
  k1 = zero(rate_prototype); k2 = zero(rate_prototype); k3 = zero(rate_prototype)

  H02 = zero(u)
  H03 = zero(u)
  H12 = Vector{typeof(u)}()
  H13 = Vector{typeof(u)}()
  H22 = Vector{typeof(u)}()
  H23 = Vector{typeof(u)}()

  for k=1:m
    push!(H12,zero(u))
    push!(H13,zero(u))
    push!(H22,zero(u))
    push!(H23,zero(u))
  end

  tmp1 = zero(rate_prototype)
  tmpg = zero(noise_rate_prototype)

  uhat = copy(uprev)
  tmp = zero(u)
  resids = zero(u)

  DRI1Cache(u,uprev,uhat,_dW,_dZ,chi1,Ihat2,tab,g1,g2,g3,k1,k2,k3,H02,H03,H12,H13,H22,H23,tmp1,tmpg,tmp,resids)
end




function alg_cache(alg::RDI4WM,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,
                   uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
    _dZ = copy(ΔW)
    chi1 = copy(ΔW)
  else
    _dW = zero(ΔW)
    _dZ = zero(ΔW)
    chi1 = zero(ΔW)
  end
  m = length(ΔW)
  Ihat2 = zeros(eltype(ΔW), m, m)
  tab = RDI4WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  g1 = zero(noise_rate_prototype)
  g2 = [zero(noise_rate_prototype) for k=1:m]
  g3 = [zero(noise_rate_prototype) for k=1:m]
  k1 = zero(rate_prototype); k2 = zero(rate_prototype); k3 = zero(rate_prototype)

  H02 = zero(u)
  H03 = zero(u)
  H12 = Vector{typeof(u)}()
  H13 = Vector{typeof(u)}()
  H22 = Vector{typeof(u)}()
  H23 = Vector{typeof(u)}()

  for k=1:m
    push!(H12,zero(u))
    push!(H13,zero(u))
    push!(H22,zero(u))
    push!(H23,zero(u))
  end

  tmp1 = zero(rate_prototype)
  tmpg = zero(noise_rate_prototype)

  uhat = copy(uprev)
  tmp = zero(u)
  resids = zero(u)

  DRI1Cache(u,uprev,uhat,_dW,_dZ,chi1,Ihat2,tab,g1,g2,g3,k1,k2,k3,H02,H03,H12,H13,H22,H23,tmp1,tmpg,tmp,resids)
end



# Roessler SRK for first order weak approx
struct RDI1WMConstantCache{T,T2} <: StochasticDiffEqConstantCache
  # hard-coded version
  a021::T

  #a121::T

  #a221::T

  b021::T

  #b121::T

  #b221::T

  α1::T
  α2::T

  c02::T2

  #c11::T2
  #c12::T2

  #c21::T2
  #c22::T2

  beta11::T
  #beta12::T

  #beta21::T
  #beta22::T

  #beta31::T
  #beta32::T

  #beta41::T
  #beta42::T

  #quantile(Normal(),1/6)
  NORMAL_ONESIX_QUANTILE::T
end



function RDI1WMConstantCache(::Type{T}, ::Type{T2}) where {T,T2}

  a021 = convert(T, 2//3)

  b021 = convert(T, 2//3)

  α1 = convert(T, 1//4)
  α2 = convert(T, 3//4)

  c02 = convert(T2, 2//3)

  beta11 = convert(T, 1)

  NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

  RDI1WMConstantCache(a021,b021,α1,α2,c02,beta11,NORMAL_ONESIX_QUANTILE)
end


function alg_cache(alg::RDI1WM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  RDI1WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


@cache struct RDI1WMCache{uType,randType,MType1,tabType,rateNoiseType,rateType,possibleRateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType

  _dW::randType
  _dZ::randType
  chi1::randType
  Ihat2::MType1

  tab::tabType

  g1::rateNoiseType

  k1::rateType
  k2::rateType

  H02::uType

  tmp1::possibleRateType
end


function alg_cache(alg::RDI1WM,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,
                   uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
    _dZ = copy(ΔW)
    chi1 = copy(ΔW)
  else
    _dW = zero(ΔW)
    _dZ = zero(ΔW)
    chi1 = zero(ΔW)
  end
  m = length(ΔW)
  Ihat2 = zeros(eltype(ΔW), m, m)
  tab = RDI1WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  g1 = zero(noise_rate_prototype)
  k1 = zero(rate_prototype); k2 = zero(rate_prototype)

  H02 = zero(u)

  tmp1 = zero(rate_prototype)

  RDI1WMCache(u,uprev,_dW,_dZ,chi1,Ihat2,tab,g1,k1,k2,H02,tmp1)
end



# Stratonovich sense

struct RSConstantCache{T,T2} <: StochasticDiffEqConstantCache
  # hard-coded version
  a021::T
  a031::T
  a032::T

  a131::T
  a141::T

  b031::T
  b032::T

  b121::T
  b131::T
  b132::T
  b141::T
  b142::T
  b143::T

  b221::T
  b231::T

  b331::T
  b332::T
  b341::T
  b342::T


  α1::T
  α2::T
  α3::T
  α4::T

  c02::T2
  c03::T2

  c13::T2
  c14::T2

  beta11::T
  beta12::T
  beta13::T
  beta14::T

  beta22::T
  beta23::T

  #quantile(Normal(),1/6)
  NORMAL_ONESIX_QUANTILE::T
end


function RS1ConstantCache(T::Type, T2::Type)
  a021 = convert(T, 0)
  a031 = convert(T, 1)
  a032 = convert(T, 0)

  a131 = convert(T, 1)
  a141 = convert(T, 1)

  b031 = convert(T, 1//4)
  b032 = convert(T, 3//4)

  b121 = convert(T, 2//3)
  b131 = convert(T, 1//12)
  b132 = convert(T, 1//4)
  b141 = convert(T, -5//4)
  b142 = convert(T, 1//4)
  b143 = convert(T, 2)

  b221 = convert(T, 1)
  b231 = convert(T, -1)

  b331 =  convert(T, 1//4)
  b332 =  convert(T, 3//4)
  b341 =  convert(T, 1//4)
  b342 =  convert(T, 3//4)

  α1 = convert(T, 0)
  α2 = convert(T, 0)
  α3 = convert(T, 1//2)
  α4 = convert(T, 1//2)

  c02 = convert(T2, 0)
  c03 = convert(T2, 1)

  c13 = convert(T2, 1)
  c14 = convert(T2, 1)

  beta11 = convert(T, 1//8)
  beta12 = convert(T, 3//8)
  beta13 = convert(T, 3//8)
  beta14 = convert(T, 1//8)

  beta22 = convert(T, -1//4)
  beta23 = convert(T, 1//4)

  NORMAL_ONESIX_QUANTILE = convert(T,-0.9674215661017014)

  RSConstantCache(a021,a031,a032,a131,a141,b031,b032,b121,b131,b132,b141,b142,b143,b221,b231,b331,b332,b341,b342,α1,α2,α3,α4,c02,c03,c13,c14,beta11,beta12,beta13,beta14,beta22,beta23,NORMAL_ONESIX_QUANTILE)
end


function alg_cache(alg::RS1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  RS1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function RS2ConstantCache(T::Type, T2::Type)
  a021 = convert(T, 2//3)
  a031 = convert(T, 1//6)
  a032 = convert(T, 1//2)

  a131 = convert(T, 1)
  a141 = convert(T, 1)

  b031 = convert(T, 1//4)
  b032 = convert(T, 3//4)

  b121 = convert(T, 2//3)
  b131 = convert(T, 1//12)
  b132 = convert(T, 1//4)
  b141 = convert(T, -5//4)
  b142 = convert(T, 1//4)
  b143 = convert(T, 2)

  b221 = convert(T, 1)
  b231 = convert(T, -1)

  b331 =  convert(T, 1//4)
  b332 =  convert(T, 3//4)
  b341 =  convert(T, 1//4)
  b342 =  convert(T, 3//4)

  α1 = convert(T, 1//4)
  α2 = convert(T, 1//4)
  α3 = convert(T, 1//2)
  α4 = convert(T, 0)

  c02 = convert(T2, 2//3)
  c03 = convert(T2, 2//3)

  c13 = convert(T2, 1)
  c14 = convert(T2, 1)

  beta11 = convert(T, 1//8)
  beta12 = convert(T, 3//8)
  beta13 = convert(T, 3//8)
  beta14 = convert(T, 1//8)

  beta22 = convert(T, -1//4)
  beta23 = convert(T, 1//4)

  NORMAL_ONESIX_QUANTILE = convert(T,-0.9674215661017014)

  RSConstantCache(a021,a031,a032,a131,a141,b031,b032,b121,b131,b132,b141,b142,b143,b221,b231,b331,b332,b341,b342,α1,α2,α3,α4,c02,c03,c13,c14,beta11,beta12,beta13,beta14,beta22,beta23,NORMAL_ONESIX_QUANTILE)
end


function alg_cache(alg::RS2,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  RS2ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


@cache struct RSCache{uType,randType,MType1,tabType,rateNoiseType,rateType,possibleRateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType

  _dW::randType
  _dZ::randType
  chi1::randType
  Ihat2::MType1

  tab::tabType

  g1::rateNoiseType
  g2::Vector{rateNoiseType}
  g3::Vector{rateNoiseType}
  g4::Vector{rateNoiseType}

  k1::rateType
  k2::rateType
  k3::rateType

  H02::uType
  H03::uType
  H12::Vector{uType}
  H13::Vector{uType}
  H14::Vector{uType}
  H22::Vector{uType}
  H23::Vector{uType}

  tmp1::possibleRateType
  tmpg::rateNoiseType

end

function alg_cache(alg::RS1,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,
                   uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
    _dZ = copy(ΔW)
    chi1 = copy(ΔW)
  else
    _dW = zero(ΔW)
    _dZ = zero(ΔW)
    chi1 = zero(ΔW)
  end
  m = length(ΔW)
  Ihat2 = zeros(eltype(ΔW), m, m)
  tab = RS1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  g1 = zero(noise_rate_prototype)
  g2 = [zero(noise_rate_prototype) for k=1:m]
  g3 = [zero(noise_rate_prototype) for k=1:m]
  g4 = [zero(noise_rate_prototype) for k=1:m]
  k1 = zero(rate_prototype); k2 = zero(rate_prototype); k3 = zero(rate_prototype)

  H02 = zero(u)
  H03 = zero(u)
  H12 = Vector{typeof(u)}()
  H13 = Vector{typeof(u)}()
  H14 = Vector{typeof(u)}()
  H22 = Vector{typeof(u)}()
  H23 = Vector{typeof(u)}()

  for k=1:m
    push!(H12,zero(u))
    push!(H13,zero(u))
    push!(H14,zero(u))
    push!(H22,zero(u))
    push!(H23,zero(u))
  end

  tmp1 = zero(rate_prototype)
  tmpg = zero(noise_rate_prototype)

  RSCache(u,uprev,_dW,_dZ,chi1,Ihat2,tab,g1,g2,g3,g4,k1,k2,k3,H02,H03,H12,H13,H14,H22,H23,tmp1,tmpg)
end


function alg_cache(alg::RS2,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,
                   uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
    _dZ = copy(ΔW)
    chi1 = copy(ΔW)
  else
    _dW = zero(ΔW)
    _dZ = zero(ΔW)
    chi1 = zero(ΔW)
  end
  m = length(ΔW)
  Ihat2 = zeros(eltype(ΔW), m, m)
  tab = RS2ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  g1 = zero(noise_rate_prototype)
  g2 = [zero(noise_rate_prototype) for k=1:m]
  g3 = [zero(noise_rate_prototype) for k=1:m]
  g4 = [zero(noise_rate_prototype) for k=1:m]
  k1 = zero(rate_prototype); k2 = zero(rate_prototype); k3 = zero(rate_prototype)

  H02 = zero(u)
  H03 = zero(u)
  H12 = Vector{typeof(u)}()
  H13 = Vector{typeof(u)}()
  H14 = Vector{typeof(u)}()
  H22 = Vector{typeof(u)}()
  H23 = Vector{typeof(u)}()

  for k=1:m
    push!(H12,zero(u))
    push!(H13,zero(u))
    push!(H14,zero(u))
    push!(H22,zero(u))
    push!(H23,zero(u))
  end

  tmp1 = zero(rate_prototype)
  tmpg = zero(noise_rate_prototype)

  RSCache(u,uprev,_dW,_dZ,chi1,Ihat2,tab,g1,g2,g3,g4,k1,k2,k3,H02,H03,H12,H13,H14,H22,H23,tmp1,tmpg)
end



# PL1WM
struct PL1WMConstantCache{T} <: StochasticDiffEqConstantCache
  #quantile(Normal(),1/6)
  NORMAL_ONESIX_QUANTILE::T
end

function PL1WMConstantCache(T::Type)
  NORMAL_ONESIX_QUANTILE = convert(T,-0.9674215661017014)

  PL1WMConstantCache(NORMAL_ONESIX_QUANTILE)
end

function alg_cache(alg::PL1WM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  PL1WMConstantCache(real(uBottomEltypeNoUnits))
end

struct PL1WMAConstantCache{T} <: StochasticDiffEqConstantCache
  #quantile(Normal(),1/6)
  NORMAL_ONESIX_QUANTILE::T
end

function PL1WMAConstantCache(T::Type)
  NORMAL_ONESIX_QUANTILE = convert(T,-0.9674215661017014)

  PL1WMAConstantCache(NORMAL_ONESIX_QUANTILE)
end

function alg_cache(alg::PL1WMA,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  PL1WMAConstantCache(real(uBottomEltypeNoUnits))
end


@cache struct PL1WMCache{uType,randType,rand2Type,MType1,tabType,rateNoiseType,rateType,possibleRateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType

  _dW::randType
  _dZ::rand2Type
  chi1::randType
  Ihat2::MType1

  tab::tabType

  g1::rateNoiseType

  k1::rateType
  k2::rateType

  Y::uType
  Yp::Vector{uType}
  Ym::Vector{uType}

  tmp1::possibleRateType
  tmpg1::rateNoiseType
  tmpg2::rateNoiseType
  Ulp::uType
  Ulm::uType
end

function alg_cache(alg::PL1WM,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,
                   uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
    _dZ = copy(ΔZ)
    chi1 = copy(ΔW)
  else
    _dW = zero(ΔW)
    _dZ = zero(ΔZ)
    chi1 = zero(ΔW)
  end
  m = length(ΔW)
  Ihat2 = zeros(eltype(ΔW), m, m)
  tab = PL1WMConstantCache(real(uBottomEltypeNoUnits))
  g1 = zero(noise_rate_prototype)

  k1 = zero(rate_prototype); k2 = zero(rate_prototype);

  Y = zero(u)
  Yp = Vector{typeof(u)}()
  Ym = Vector{typeof(u)}()

  for k=1:m
    push!(Yp,zero(u))
    push!(Ym,zero(u))
  end

  tmp1 = zero(rate_prototype)
  tmpg1 = zero(noise_rate_prototype)
  tmpg2 = zero(noise_rate_prototype)

  Ulp = zero(u)
  Ulm = zero(u)

  PL1WMCache(u,uprev,_dW,_dZ,chi1,Ihat2,tab,g1,k1,k2,Y,Yp,Ym,tmp1,tmpg1,tmpg2,Ulp,Ulm)
end


# additive noise
@cache struct PL1WMACache{uType,randType,tabType,rateNoiseType,rateType,possibleRateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType

  _dW::randType
  chi1::randType
  tab::tabType

  g1::rateNoiseType
  k1::rateType
  k2::rateType

  Y::uType
  tmp1::possibleRateType
end


function alg_cache(alg::PL1WMA,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,
                   uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
    chi1 = copy(ΔW)
  else
    _dW = zero(ΔW)
    chi1 = zero(ΔW)
  end

  tab = PL1WMConstantCache(real(uBottomEltypeNoUnits))
  g1 = zero(noise_rate_prototype)

  k1 = zero(rate_prototype); k2 = zero(rate_prototype);

  Y = zero(u)

  tmp1 = zero(rate_prototype)

  PL1WMACache(u,uprev,_dW,chi1,tab,g1,k1,k2,Y,tmp1)
end
