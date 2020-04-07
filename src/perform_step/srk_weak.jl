
@muladd function perform_step!(integrator,cache::DRI1ConstantCache,f=integrator.f)
  @unpack a021,a031,a032,a121,a131,b021,b031,b121,b131,b221,b222,b223,b231,b232,b233,α1,α2,α3,c02,c03,c12,c13,beta11,beta12,beta13,beta22,beta23,beta31,beta32,beta33,beta42,beta43,NORMAL_ONESIX_QUANTILE = cache
  @unpack t,dt,uprev,u,W,p = integrator

  # define three-point distributed random variables
  dW_scaled = W.dW / sqrt(dt)
  _dW = map(x -> calc_threepoint_random(integrator, NORMAL_ONESIX_QUANTILE, x), dW_scaled)
  Ihat2_diag = map(x -> (x^2-dt)/2, _dW)
  if !(typeof(W.dW) <: Number)
    m = length(W.dW)
    # define two-point distributed random variables
    _dZ = map(x -> calc_twopoint_random(integrator, x),  W.dZ)
    Ihat2 = zeros(eltype(W.dZ), m, m) # I^_(k,l)
    for k = 1:m
      for l = 1:m
        if k<l
          Ihat2[k, l] = (_dW[k]*_dW[l]-integrator.sqdt*_dZ[k])/2
        elseif l<k
          Ihat2[k, l] = (_dW[k]*_dW[l]+integrator.sqdt*_dZ[l])/2
        else k==l
          Ihat2[k, k] = Ihat2_diag[k]
        end
      end
    end
  end

  # compute stage values
  k1 = integrator.f(uprev,p,t)
  g1 = integrator.g(uprev,p,t)

  # H_i^(0), stage 1
  # H01 = uprev
  # H_i^(0), stage 2
  if !is_diagonal_noise(integrator.sol.prob) || typeof(W.dW) <: Number
    H02 = uprev + a021*k1*dt + b021*g1*_dW
  else
    H02 = uprev + a021*k1*dt + b021*g1.*_dW
  end

  # H_i^(0), stage 3 (requires H_i^(k) stage 2 in general)
  k2 = integrator.f(H02,p,t+c02*dt)
  H03 = uprev + a032*k2*dt + a031*k1*dt
  if !is_diagonal_noise(integrator.sol.prob) || typeof(W.dW) <: Number
    H03 += b031*g1*_dW
  else
    H03 += b031*g1.*_dW
  end

  k3 = integrator.f(H03,p,t+c03*dt)

  # H_i^(k), stage 1
  # H11 = uprev
  # H_i^(k), stage 2
  if typeof(W.dW) <: Number
    H12 = uprev + a121*k1*dt + b121*g1*integrator.sqdt
  elseif is_diagonal_noise(integrator.sol.prob)
    H12 = Vector{typeof(uprev)}[uprev .+ a121*k1*dt .+ b121*g1[k]*integrator.sqdt for k=1:m]
  else
    H12 = Vector{typeof(uprev)}[uprev .+ a121*k1*dt .+ b121*g1[:,k]*integrator.sqdt for k=1:m]
  end

  # # H_i^(k), stage 3
  if typeof(W.dW) <: Number
    H13 = uprev + a131*k1*dt + b131*g1*integrator.sqdt
  elseif is_diagonal_noise(integrator.sol.prob)
    H13 = Vector{typeof(uprev)}[uprev .+ a131*k1*dt .+ b131*g1[k]*integrator.sqdt for k=1:m]
  else
    H13 = Vector{typeof(uprev)}[uprev .+ a131*k1*dt .+ b131*g1[:,k]*integrator.sqdt for k=1:m]
  end

  # H^_i^(k), stage 1
  # H21 = uprev
  # H^_i^(k), stage 2

  if typeof(W.dW) <: Number
    g2 = integrator.g(H12,p,t+c12*dt)
    g3 = integrator.g(H13,p,t+c13*dt)
    # for m=1:  H22 = uprev
  else
    g2 = [integrator.g(H12[k],p,t+c12*dt) for k=1:m]
    g3 = [integrator.g(H13[k],p,t+c13*dt) for k=1:m]
    H22 = [uprev for k=1:m]
    # add inbounds for speed if working properly
    for k=1:m
      for l=1:m
        if k == l
          continue
        end
        if is_diagonal_noise(integrator.sol.prob)
          @.. H22[k] += (b221*g1[l]+b222*g2[l][l]+b223*g3[l][l])*Ihat2[k,l]/integrator.sqdt
        else
          H22[k] += (b221*g1[:,l]+b222*g2[l][:,l]+b223*g3[l][:,l])*Ihat2[k,l]/integrator.sqdt
        end
      end
    end
  end

  # H^_i^(k), stage 3
  # for m=1: H23 = uprev
  if !(typeof(W.dW) <: Number)
    H23 = [uprev for k=1:m]
    # add inbounds for speed if working properly
    for k=1:m
      for l=1:m
        if k == l
          continue
        end
        if is_diagonal_noise(integrator.sol.prob)
          @.. H23[k] += (b231*g1[l]+b232*g2[l][l]+b233*g3[l][l])*Ihat2[k,l]/integrator.sqdt
        else
          H23[k] += (b231*g1[:,l]+b232*g2[l][:,l]+b233*g3[l][:,l])*Ihat2[k,l]/integrator.sqdt
        end
      end
    end
  end

  # add stages together Eq. (3)
  u = uprev + α1*k1*dt + α2*k2*dt + α3*k3*dt

  # add noise
  if typeof(W.dW) <: Number
    # lines 2 and 3
    u += g1*(_dW*beta11) + g2*(_dW*beta12+Ihat2_diag*beta22/integrator.sqdt) + g3*(_dW*beta13+Ihat2_diag*beta23/integrator.sqdt)
    # lines 4 and 5 are zero by construction
    # u += g1*(_dW*(beta31+beta32+beta33)+Ihat2_diag*integrator.sqdt*(beta42+beta43))

  else
    for k=1:m
      if is_diagonal_noise(integrator.sol.prob)
        u += g1.*_dW*beta11 + g2[k].*_dW*beta12 + g3[k].*_dW*beta13
        u += g2[k].*Ihat2_diag/integrator.sqdt*beta22 + g3[k].*Ihat2_diag/integrator.sqdt*beta23
        u += g1.*_dW*beta31 + integrator.g(H22[k],p,t).*(_dW*beta32 .+ integrator.sqdt*beta42) + integrator.g(H23[k],p,t).*(_dW*beta33 .+ integrator.sqdt*beta43)
      else
        # non-diag noise
        u += g1*_dW*beta11 + g2[k]*_dW*beta12 + g3[k]*_dW*beta13
        u += g2[k]*Ihat2_diag/integrator.sqdt*beta22 + g3[k]*Ihat2_diag/integrator.sqdt*beta23
        u += g1*_dW*beta31 + integrator.g(H22[k],p,t)*(_dW*beta32 .+ integrator.sqdt*beta42) + integrator.g(H23[k],p,t)*(_dW*beta33 .+ integrator.sqdt*beta43)
      end
    end
  end

  # Test EM scheme:
  # integrator.u = uprev + 1*k1*dt + 1*g1*_dW

  integrator.u = u

end
