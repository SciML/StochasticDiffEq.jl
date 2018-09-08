struct LinearInterpolationData{uType,tType} <: AbstractDiffEqInterpolation
  timeseries::uType
  ts::tType
end

DiffEqBase.interp_summary(::LinearInterpolationData) = "1st order linear"
(interp::LinearInterpolationData)(tvals,idxs,deriv,p,continuity::Symbol=:left) = sde_interpolation(tvals,interp,idxs,deriv,p,continuity)
(interp::LinearInterpolationData)(val,tvals,idxs,deriv,p,continuity::Symbol=:left) = sde_interpolation!(val,tvals,interp,idxs,deriv,p,continuity)
