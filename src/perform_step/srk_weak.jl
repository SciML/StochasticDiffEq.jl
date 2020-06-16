
@muladd function perform_step!(integrator,cache::DRI1ConstantCache,f=integrator.f)
  @unpack a021,a031,a032,a121,a131,b021,b031,b121,b131,b221,b222,b223,b231,b232,b233,α1,α2,α3,c02,c03,c12,c13,beta11,beta12,beta13,beta22,beta23,beta31,beta32,beta33,beta42,beta43,NORMAL_ONESIX_QUANTILE = cache
  @unpack t,dt,uprev,u,W,p = integrator

  # define three-point distributed random variables
  dW_scaled = W.dW / sqrt(dt)
  _dW = map(x -> calc_threepoint_random(integrator, NORMAL_ONESIX_QUANTILE, x), dW_scaled)
  chi1 = map(x -> (x^2-dt)/2, _dW) # diagonal of Ihat2
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
          Ihat2[k, k] = chi1[k]
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
    u += g1*(_dW*beta11) + g2*(_dW*beta12+chi1*beta22/integrator.sqdt) + g3*(_dW*beta13+chi1*beta23/integrator.sqdt)
    # lines 4 and 5 are zero by construction
    # u += g1*(_dW*(beta31+beta32+beta33)+chi1*integrator.sqdt*(beta42+beta43))
  else
      if is_diagonal_noise(integrator.sol.prob)
        u += g1.*_dW*beta11
        for k=1:m
          u[k] += g2[k][k]*_dW[k]*beta12+g3[k][k]*_dW[k]*beta13
          u[k] += g2[k][k]*chi1[k]*beta22/integrator.sqdt + g3[k][k]*chi1[k]*beta23/integrator.sqdt
          tmpg = integrator.g(H22[k],p,t)
          u[k] = u[k] + g1[k]*_dW[k]*beta31 + tmpg[k]*(_dW[k]*beta32 + integrator.sqdt*beta42)
          tmpg = integrator.g(H23[k],p,t)
          u[k] = u[k] + tmpg[k]*(_dW[k]*beta33 + integrator.sqdt*beta43)
        end
      else
        # non-diag noise
        for k=1:m
          g1k = @view g1[:,k]
          g2k = @view g2[k][:,k]
          g3k = @view g3[k][:,k]
          @.. u = u + g1k*_dW[k]*(beta11+beta31)
          @.. u = u + g2k*_dW[k]*beta12 + g3k*_dW[k]*beta13
          @.. u = u + g2k*chi1[k]*beta22/integrator.sqdt + g3k*chi1[k]*beta23/integrator.sqdt
          tmpg = integrator.g(H22[k],p,t)
          @.. u = u + tmpgk*_dW[k]*beta32 + tmpgk*integrator.sqdt*beta42
          tmpg = integrator.g(H23[k],p,t)
          @.. u = u + tmpgk*_dW[k]*beta33 + tmpgk*integrator.sqdt*beta43
        end
      end
  end

  if integrator.opts.adaptive && (typeof(W.dW) <: Number || is_diagonal_noise(integrator.sol.prob) || m==1)

    # schemes with lower convergence order
    if c03!=0.0
      # scheme has det. conv. order 3 and we look for det. conv. order 2 scheme
      rat = c02/c03
      δ₁ = integrator.opts.delta*(rat-1)
      δ₂ = -integrator.opts.delta*rat
      uhat = uprev + ((α1+δ₁)*k1 + (α2+integrator.opts.delta)*k2 + (α3+δ₂)*k3)*dt
    else
      # check against EM
      uhat = uprev + k1*dt
    end
    if is_diagonal_noise(integrator.sol.prob)
      uhat += g1 .* _dW
    else
      uhat += g1 * _dW
    end

    Eprev = Statistics.mean(uprev,dims=2)
    E₁ = Statistics.mean(u,dims=2)
    E₂ = Statistics.mean(uhat,dims=2)

    resids = calculate_residuals(E₁-E₂,Eprev,E₁,integrator.opts.abstol,integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(resids, t)
  end

  integrator.u = u

end



@muladd function perform_step!(integrator,cache::DRI1Cache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  @unpack _dW,_dZ,chi1,Ihat2,tab,g1,g2,g3,k1,k2,k3,H02,H03,H12,H13,H22,H23,tmp1,tmpg,uhat,Eprev,E₁,E₂,resids = cache
  @unpack a021,a031,a032,a121,a131,b021,b031,b121,b131,b221,b222,b223,b231,b232,b233,α1,α2,α3,c02,c03,c12,c13,beta11,beta12,beta13,beta22,beta23,beta31,beta32,beta33,beta42,beta43,NORMAL_ONESIX_QUANTILE = cache.tab

  m = length(W.dW)

  if typeof(W.dW) <: Union{SArray,Number}
    # tbd
    _dW = map(x -> calc_threepoint_random(integrator, NORMAL_ONESIX_QUANTILE, x), W.dW / sqrt(dt))
  else
    # define three-point distributed random variables
    @.. chi1 = W.dW / sqrt(dt)
    calc_threepoint_random!(_dW, integrator, NORMAL_ONESIX_QUANTILE, chi1)
    map!(x -> (x^2-dt)/2, chi1, _dW)

    # define two-point distributed random variables
    calc_twopoint_random!(_dZ,integrator,W.dZ)

    for k = 1:m
      for l = 1:m
        if k<l
          Ihat2[k, l] = (_dW[k]*_dW[l]-integrator.sqdt*_dZ[k])/2
        elseif l<k
          Ihat2[k, l] = (_dW[k]*_dW[l]+integrator.sqdt*_dZ[l])/2
        else k==l
          Ihat2[k, k] = chi1[k]
        end
      end
    end
  end

  # compute stage values
  integrator.f(k1,uprev,p,t)
  integrator.g(g1,uprev,p,t)

  # H_i^(0), stage 1
  # H01 = uprev
  # H_i^(0), stage 2
  if is_diagonal_noise(integrator.sol.prob)
    @.. H02 = uprev + a021*k1*dt + b021*g1*_dW
  else
    mul!(tmp1,g1,_dW)
    @.. H02 = uprev + dt*a021*k1 + b021*tmp1
  end

  integrator.f(k2,H02,p,t+c02*dt)
  if is_diagonal_noise(integrator.sol.prob)
    @.. H03 = uprev + a032*k2*dt + a031*k1*dt + b031*g1*_dW
  else
    @.. H03 = uprev + a032*k2*dt + a031*k1*dt + b031*tmp1
  end

  integrator.f(k3,H03,p,t+c03*dt)

  # H_i^(k), stages
  # H11 = uprev
  for k=1:m
    if is_diagonal_noise(integrator.sol.prob)
      @.. H12[k] = uprev + a121*k1*dt + b121*g1[k]*integrator.sqdt
      @.. H13[k] = uprev + a131*k1*dt + b131*g1[k]*integrator.sqdt
    else
      g1k = @view g1[:,k]
      @.. H12[k] = uprev + a121*k1*dt + b121*g1k*integrator.sqdt
      @.. H13[k] = uprev + a131*k1*dt + b131*g1k*integrator.sqdt
    end
  end
  # H^_i^(k), stages

  for k=1:m
    integrator.g(g2[k],H12[k],p,t+c12*dt)
    integrator.g(g3[k],H13[k],p,t+c13*dt)

    H22[k] = uprev
    H23[k] = uprev
  end

  for k=1:m
    for l=1:m
      if k == l
        continue
      end
      if is_diagonal_noise(integrator.sol.prob)
        @.. H22[k] = H22[k]+(b221*g1[l]+b222*g2[l][l]+b223*g3[l][l])*Ihat2[k,l]/integrator.sqdt
        @.. H23[k] = H23[k]+(b231*g1[l]+b232*g2[l][l]+b233*g3[l][l])*Ihat2[k,l]/integrator.sqdt
      else
        g1l = @view g1[:,l]
        g2l = @view g2[l][:,l]
        g3l = @view g3[l][:,l]
        @.. H22[k] = H22[k]+(b221*g1l+b222*g2l+b223*g3l)*Ihat2[k,l]/integrator.sqdt
        @.. H23[k] = H23[k]+(b231*g1l+b232*g2l+b233*g3l)*Ihat2[k,l]/integrator.sqdt
      end
    end
  end

  # add stages together Eq. (3)
  @.. u = uprev + α1*k1*dt + α2*k2*dt + α3*k3*dt

  # add noise
  if is_diagonal_noise(integrator.sol.prob)
      @.. u = u + g1*_dW*beta11
      for k=1:m
        u[k] = u[k] + g2[k][k]*_dW[k]*beta12 + g3[k][k]*_dW[k]*beta13
        u[k] = u[k]  + g2[k][k]*chi1[k]*beta22/integrator.sqdt + g3[k][k]*chi1[k]*beta23/integrator.sqdt
        if m > 1
          integrator.g(tmpg,H22[k],p,t)
          u[k] = u[k] + g1[k]*_dW[k]*beta31 + tmpg[k]*(_dW[k]*beta32 + integrator.sqdt*beta42)
          integrator.g(tmpg,H23[k],p,t)
          u[k] = u[k] + tmpg[k]*(_dW[k]*beta33 + integrator.sqdt*beta43)
        end
      end
  else
    for k=1:m
      g1k = @view g1[:,k]
      g2k = @view g2[k][:,k]
      g3k = @view g3[k][:,k]
      tmpgk = @view tmpg[:,k]
      @.. u = u + g1k*_dW[k]*(beta11+beta31)
      @.. u = u + g2k*_dW[k]*beta12 + g3k*_dW[k]*beta13
      @.. u = u + g2k*chi1[k]*beta22/integrator.sqdt + g3k*chi1[k]*beta23/integrator.sqdt
      integrator.g(tmpg,H22[k],p,t)
      @.. u = u + tmpgk*_dW[k]*beta32 + tmpgk*integrator.sqdt*beta42
      integrator.g(tmpg,H23[k],p,t)
      @.. u = u + tmpgk*_dW[k]*beta33 + tmpgk*integrator.sqdt*beta43
    end
  end

  if integrator.opts.adaptive && (typeof(W.dW) <: Number || is_diagonal_noise(integrator.sol.prob) || m==1)

    # schemes with lower convergence order
    if c03!=0.0
      # scheme has det. conv. order 3 and we look for det. conv. order 2 scheme
      rat = c02/c03
      δ₁ = integrator.opts.delta*(rat-1)
      δ₂ = -integrator.opts.delta*rat
      @.. uhat = uprev + ((α1+δ₁)*k1 + (α2+integrator.opts.delta)*k2 + (α3+δ₂)*k3)*dt
    else
      # check against EM
      @.. uhat = uprev + k1*dt
    end
    if is_diagonal_noise(integrator.sol.prob)
      @.. uhat = uhat + g1 .* _dW
    else
      @.. uhat = uhat + g1 * _dW
    end

    Statistics.mean!(Eprev,uprev)
    Statistics.mean!(E₁,u)
    Statistics.mean!(E₂,uhat)

    @.. E₂ = E₁-E₂

    calculate_residuals!(resids, E₂, Eprev, E₁, integrator.opts.abstol,
                                integrator.opts.reltol,integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(resids, t)
  end

end



# Roessler SRK for first order weak approx
@muladd function perform_step!(integrator,cache::RDI1WMConstantCache,f=integrator.f)
  @unpack a021,b021,α1,α2,c02,beta11,NORMAL_ONESIX_QUANTILE = cache
  @unpack t,dt,uprev,u,W,p = integrator

  # define three-point distributed random variables
  dW_scaled = W.dW / sqrt(dt)
  _dW = map(x -> calc_threepoint_random(integrator, NORMAL_ONESIX_QUANTILE, x), dW_scaled)
  chi1 = map(x -> (x^2-dt)/2, _dW) # diagonal of Ihat2
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
          Ihat2[k, k] = chi1[k]
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
  k2 = integrator.f(H02,p,t+c02*dt)

  # add stages together Eq. (3)
  u = uprev + α1*k1*dt + α2*k2*dt

  # add noise
  if typeof(W.dW) <: Number
    # lines 2 and 3
    u += g1*(_dW*beta11)
  else
      if is_diagonal_noise(integrator.sol.prob)
        u += g1.*_dW*beta11
      else
        # non-diag noise
        for k=1:m
          g1k = @view g1[:,k]
          @.. u = u + g1k*_dW[k]*(beta11)
        end
      end
  end

  integrator.u = u

end



@muladd function perform_step!(integrator,cache::RDI1WMCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  @unpack _dW,_dZ,chi1,Ihat2,tab,g1,k1,k2,H02,tmp1 = cache
  @unpack a021,b021,α1,α2,c02,beta11,NORMAL_ONESIX_QUANTILE = cache.tab

  m = length(W.dW)

  if typeof(W.dW) <: Union{SArray,Number}
    # tbd
    _dW = map(x -> calc_threepoint_random(integrator, NORMAL_ONESIX_QUANTILE, x), W.dW / sqrt(dt))
  else
    # define three-point distributed random variables
    @.. chi1 = W.dW / sqrt(dt)
    calc_threepoint_random!(_dW, integrator, NORMAL_ONESIX_QUANTILE, chi1)
    map!(x -> (x^2-dt)/2, chi1, _dW)

    # define two-point distributed random variables
    calc_twopoint_random!(_dZ,integrator,W.dZ)

    for k = 1:m
      for l = 1:m
        if k<l
          Ihat2[k, l] = (_dW[k]*_dW[l]-integrator.sqdt*_dZ[k])/2
        elseif l<k
          Ihat2[k, l] = (_dW[k]*_dW[l]+integrator.sqdt*_dZ[l])/2
        else k==l
          Ihat2[k, k] = chi1[k]
        end
      end
    end
  end

  # compute stage values
  integrator.f(k1,uprev,p,t)
  integrator.g(g1,uprev,p,t)

  # H_i^(0), stage 1
  # H01 = uprev
  # H_i^(0), stage 2
  if is_diagonal_noise(integrator.sol.prob)
    @.. H02 = uprev + a021*k1*dt + b021*g1*_dW
  else
    mul!(tmp1,g1,_dW)
    @.. H02 = uprev + dt*a021*k1 + b021*tmp1
  end

  integrator.f(k2,H02,p,t+c02*dt)

  # add stages together Eq. (3)
  @.. u = uprev + α1*k1*dt + α2*k2*dt

  # add noise
  if is_diagonal_noise(integrator.sol.prob)
      @.. u = u + g1*_dW*beta11
  else
    for k=1:m
      g1k = @view g1[:,k]
      @.. u = u + g1k*_dW[k]*beta11
    end
  end
end
