immutable LinearInterpolationData{uType,tType} <: Function
  timeseries::uType
  ts::tType
end

(interp::LinearInterpolationData)(tvals) = sde_interpolation(tvals,interp)
