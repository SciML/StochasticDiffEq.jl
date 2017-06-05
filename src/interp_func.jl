immutable LinearInterpolationData{uType,tType} <: AbstractDiffEqInterpolation
  timeseries::uType
  ts::tType
end

DiffEqBase.interp_summary(::LinearInterpolationData) = "First Order Linear"
(interp::LinearInterpolationData)(tvals,idxs,deriv) = sde_interpolation(tvals,interp,idxs,deriv)
(interp::LinearInterpolationData)(val,tvals,idxs,deriv) = sde_interpolation!(val,tvals,interp,idxs,deriv)
