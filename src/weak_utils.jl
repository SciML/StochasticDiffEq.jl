@muladd function calc_twopoint_random!(_dW, sqdt, dW)
  @.. _dW = ifelse(sign(dW) > false, sqdt, -sqdt)
  return nothing
end


@muladd function calc_twopoint_random(sqdt, dW)
  return sign(dW) > 0.0 ? sqdt : -sqdt
end


@muladd function calc_threepoint_random!(_dW, integrator, quantile, dW_scaled)
  for (index, value) in enumerate(dW_scaled)
    if value < quantile
      _dW[index] = -sqrt(3*integrator.dt)
    elseif value > -quantile
      _dW[index] = sqrt(3*integrator.dt)
    else
      _dW[index] = zero(integrator.dt)
    end
  end
end


@muladd function calc_threepoint_random(integrator, quantile, dW_scaled)
  if dW_scaled < quantile
    _dW = -sqrt(3*integrator.dt)
  elseif dW_scaled > -quantile
    _dW = sqrt(3*integrator.dt)
  else
    _dW = zero(integrator.dt)
  end
  _dW
end
