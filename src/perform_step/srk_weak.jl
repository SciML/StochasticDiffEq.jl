

@muladd function perform_step!(integrator,cache::DRI1ConstantCache,f=integrator.f)
  @unpack a021,a031,a032,a121,a131,b021,b031,b121,b131,b221,b222,b223,b231,b232,b233,α1,α2,α3,c02,c03,c12,c13,beta11,beta12,beta13,beta22,beta23,beta31,beta32,beta33,beta42,beta43 = cache
  @unpack t,dt,uprev,u,W,p = integrator

  # define random processes
  # table method for sampling
  sqrt3dt = sqrt(3one(eltype(W.dW)))*integrator.sqdt
  I1_outcomes = [-sqrt3dt, sqrt3dt, zero(eltype(W.dW)), zero(eltype(W.dW)), zero(eltype(W.dW)), zero(eltype(W.dW))]
  if typeof(W.dW) <: Number
    m = 1  # dim of Wiener process
    picks = rand(1:6) # to have it a simple float instead of an array
  else
    m = length(W.dW)
    picks = rand(1:6, m)
  end
  # three-point distributed random variables
  Ihat1 = I1_outcomes[picks] # I^_(k)

  # two-point distributed random variables
  if typeof(W.dW) <: Number
    Ihat2 = (Ihat1^2-dt)/2
  else # m>1
    diagIhat2 = (Ihat1.^2 .-dt)/2
    picks = rand(1:2, m-1)
    I2_outcomes = [-integrator.sqdt, integrator.sqdt]
    Itilde = I2_outcomes[picks] # I~_(k)
    Ihat2 = zeros(eltype(W.dW), m, m) # I^_(k,l)
    for k = 1:m
      for l = 1:m
        if l<k
          Ihat2[k, l] = (Ihat1[k]*Ihat1[l]+integrator.sqdt*Itilde[l])/2
        elseif k<l
          Ihat2[k, l] = (Ihat1[k]*Ihat1[l]-integrator.sqdt*Itilde[k])/2
        else k==l
          Ihat2[k, k] = diagIhat2[k]
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
  # @show !is_diagonal_noise(integrator.sol.prob), typeof(W.dW)<: Number
  if !is_diagonal_noise(integrator.sol.prob) || typeof(W.dW) <: Number
    H02 = uprev + a021*k1*dt + b021*g1*Ihat1
  else
    H02 = uprev + a021*k1*dt + b021*g1.*Ihat1
  end

  # H_i^(0), stage 3 (requires H_i^(k) stage 2 in general)
  k2 = integrator.f(H02,p,t+c02*dt)
  H03 = uprev + a032*k2*dt + a031*k1*dt
  if !is_diagonal_noise(integrator.sol.prob) || typeof(W.dW) <: Number
    H03 += b031*g1*Ihat1
  else
    H03 += b031*g1.*Ihat1
  end

  k3 = integrator.f(H03,p,t+c03*dt)

  # H_i^(k), stage 1
  # H11 = uprev
  # H_i^(k), stage 2
  if typeof(W.dW) <: Number
    H12 = uprev + a121*k1*dt + b121*g1*integrator.sqdt
  elseif is_diagonal_noise(integrator.sol.prob)
    H12 = [uprev + a121*k1*dt .+ b121*g1[k]*integrator.sqdt for k=1:m]
  else
    H12 = [uprev + a121*k1*dt + b121*g1[:,k]*integrator.sqdt for k=1:m]
  end

  # # H_i^(k), stage 3
  if typeof(W.dW) <: Number
    H13 = uprev + a131*k1*dt + b131*g1*integrator.sqdt
  elseif is_diagonal_noise(integrator.sol.prob)
    H13 = [uprev + a131*k1*dt .+ b131*g1[k]*integrator.sqdt for k=1:m]
  else
    H13 = [uprev + a131*k1*dt + b131*g1[:,k]*integrator.sqdt for k=1:m]
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
    u += g1*(Ihat1*beta11) + g2*(Ihat1*beta12+Ihat2*beta22/integrator.sqdt) + g3*(Ihat1*beta13+Ihat2*beta23/integrator.sqdt)
    # lines 4 and 5 are zero by construction
    # u += g1*(Ihat1*(beta31+beta32+beta33)+Ihat2*integrator.sqdt*(beta42+beta43))

  else

    for k=1:m
      if is_diagonal_noise(integrator.sol.prob)
        u += g1.*Ihat1*beta11 + g2[k].*Ihat1*beta12 + g3[k].*Ihat1*beta13
        u += g2[k].*diagIhat2/integrator.sqdt*beta22 + g3[k].*diagIhat2/integrator.sqdt*beta23
        u += g1.*Ihat1*beta31 + integrator.g(H22[k],p,t).*(Ihat1*beta32 .+ integrator.sqdt*beta42) + integrator.g(H23[k],p,t).*(Ihat1*beta33 .+ integrator.sqdt*beta43)     
      else
        # non-diag noise
        u += g1*Ihat1*beta11 + g2[k]*Ihat1*beta12 + g3[k]*Ihat1*beta13
        u += g2[k]*diagIhat2/integrator.sqdt*beta22 + g3[k]*diagIhat2/integrator.sqdt*beta23
        u += g1*Ihat1*beta31 + integrator.g(H22[k],p,t)*(Ihat1*beta32 .+ integrator.sqdt*beta42) + integrator.g(H23[k],p,t)*(Ihat1*beta33 .+ integrator.sqdt*beta43)
      end
    end
  end

  # Test EM scheme:
  #integrator.u = uprev + 1*k1*dt + 1*g1*Ihat1

  integrator.u = u
end
