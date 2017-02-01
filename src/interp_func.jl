immutable LinearInterpolationData{uType,tType} <: Function
  timeseries::uType
  ts::tType
end

(interp::LinearInterpolationData)(tvals,idxs,deriv) = sde_interpolation(tvals,interp,idxs,deriv)
(interp::LinearInterpolationData)(val,tvals,idxs,deriv) = sde_interpolation!(val,tvals,interp,idxs,deriv)
