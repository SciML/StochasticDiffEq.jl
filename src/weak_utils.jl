@muladd function calc_twopoint_random!(_dW, integrator, dW)
    @.. _dW = ifelse(sign(dW) > 0.0, integrator.sqdt, -integrator.sqdt)
end

@muladd function calc_twopoint_random(integrator, dW)
  return sign(dW) > 0.0 ? integrator.sqdt : -integrator.sqdt
end
