mutable struct TimeGradientWrapper{VFType,uType} <: Function
  vf::VFType
  uprev::uType
end
(p::TimeGradientWrapper)(t) = (du2 = similar(p.uprev); p.vf(t,p.uprev,du2); du2)
(p::TimeGradientWrapper)(du2,t) = p.vf(t,p.uprev,du2)

mutable struct UJacobianWrapper{VFRType,tType} <: Function
  vfr::VFRType
  t::tType
end

(p::UJacobianWrapper)(du1,uprev) = p.vfr(p.t,uprev,du1)
(p::UJacobianWrapper)(uprev) = (du1 = similar(uprev); p.vfr(p.t,uprev,du1); du1)

mutable struct TimeDerivativeWrapper{F,uType} <: Function
  f::F
  u::uType
end
(p::TimeDerivativeWrapper)(t) = p.f(t,p.u)

mutable struct UDerivativeWrapper{F,tType} <: Function
  f::F
  t::tType
end
(p::UDerivativeWrapper)(u) = p.f(p.t,u)
