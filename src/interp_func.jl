struct LinearInterpolationData{uType, tType} <: AbstractDiffEqInterpolation
    timeseries::uType
    ts::tType
end

DiffEqBase.interp_summary(::LinearInterpolationData) = "1st order linear"
function (interp::LinearInterpolationData)(tvals, idxs, deriv, p,
                                           continuity::Symbol = :left)
    sde_interpolation(tvals, interp, idxs, deriv, p, continuity)
end
function (interp::LinearInterpolationData)(val, tvals, idxs, deriv, p,
                                           continuity::Symbol = :left)
    sde_interpolation!(val, tvals, interp, idxs, deriv, p, continuity)
end
