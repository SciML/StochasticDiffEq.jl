### SDE Examples

f = (t,u) -> 1.01*u
σ = (t,u) -> 0.87*u
analytic = (t,u₀,W) -> u₀.*exp.(0.63155*t+0.87*W)
"""
```math
du_t = βudt + αudW_t
```
where β=1.01, α=0.87, and initial condtion u₀=1/2, with solution

```math
u(t,u₀,W_t)=u₀\\exp((α-\\frac{β^2}{2})t+βW_t)
```

"""
prob_sde_linear = SDEProblem(f,σ,1/2,analytic=analytic)

f = (t,u,du) -> begin
  for i = 1:length(u)
    du[i] = 1.01*u[i]
  end
end
σ = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = .87*u[i]
  end
end
"""
8 linear SDEs (as a 4x2 matrix):

```math
du_t = βudt + αudW_t
```
where β=1.01, α=0.87, and initial condtion u₀=1/2 with solution

```math
u(t,u₀,W_t)=u₀\\exp((α-\\frac{β^2}{2})t+βW_t)
```
"""
prob_sde_2Dlinear = SDEProblem(f,σ,ones(4,2)/2,analytic=analytic)


f = (t,u) -> -.25*u*(1-u^2)
σ = (t,u) -> .5*(1-u^2)
analytic = (t,u₀,W) -> ((1+u₀).*exp.(W)+u₀-1)./((1+u₀).*exp.(W)+1-u₀)
"""
```math
du_t = \\frac{1}{4}u(1-u^2)dt + \\frac{1}{2}(1-u^2)dW_t
```

and initial condtion u₀=1/2, with solution

```math
u(t,u₀,W_t)=\\frac{(1+u₀)\\exp(W_t)+u₀-1}{(1+u₀)\\exp(W_t)+1-u₀}
```
"""
prob_sde_cubic = SDEProblem(f,σ,1/2,analytic=analytic)

f = (t,u) -> -0.01*sin.(u).*cos.(u).^3
σ = (t,u) -> 0.1*cos.(u).^2
analytic = (t,u₀,W) -> atan.(0.1*W + tan.(u₀))
"""
```math
du_t = -\\frac{1}{100}\sin(u)\cos^3(u)dt + \\frac{1}{10}\cos^{2}(u_t) dW_t
```

and initial condition `u₀=1.0` with solution

```math
u(t,u₀,W_t)=\\arctan(\\frac{W_t}{10} + \\tan(u₀))
```
"""
prob_sde_wave = SDEProblem(f,σ,1.,analytic=analytic)

const sde_wave_α = 0.1
const sde_wave_β = 0.05
f = (t,u) -> sde_wave_β./sqrt.(1+t) - u./(2*(1+t))
σ = (t,u) -> sde_wave_α*sde_wave_β./sqrt.(1+t)
analytic = (t,u₀,W) -> u₀./sqrt.(1+t) + sde_wave_β*(t+sde_wave_α*W)./sqrt.(1+t)

"""
Additive noise problem

```math
u_t = (\\frac{β}{\\sqrt{1+t}}-\\frac{1}{2(1+t)}u_t)dt + \\frac{αβ}{\\sqrt{1+t}}dW_t
```

and initial condition u₀=1.0 with α=0.1 and β=0.05, with solution

```math
u(t,u₀,W_t)=\\frac{u₀}{\\sqrt{1+t}} + \\frac{β(t+αW_t)}{\\sqrt{1+t}}
```
"""
prob_sde_additive = SDEProblem(f,σ,1.,analytic=analytic)

const sde_wave_αvec = [0.1;0.1;0.1;0.1]
const sde_wave_βvec = [0.5;0.25;0.125;0.1115]
f = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = sde_wave_βvec[i]/sqrt(1+t) - u[i]/(2*(1+t))
  end
end

σ = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = sde_wave_αvec[i]*sde_wave_βvec[i]/sqrt(1+t)
  end
end
analytic = (t,u₀,W) -> u₀./sqrt(1+t) + sde_wave_βvec.*(t+sde_wave_αvec.*W)./sqrt(1+t)

"""
A multiple dimension extension of `additiveSDEExample`

"""
prob_sde_additivesystem = SDEProblem(f,σ,[1.;1.;1.;1.],analytic=analytic)

f = @ode_def_nohes Lorenz begin
  dx = σ*(y-x)
  dy = x*(ρ-z) - y
  dz = x*y - β*z
end σ=>10. ρ=>28. β=>2.66

σ = (t,u) -> 3.0 #Additive
"""
Lorenz Attractor with additive noise

```math
\\begin{align}
dx &= σ*(y-x)dt + αdW_t \\\\
dy &= (x*(ρ-z) - y)dt + αdW_t \\\\
dz &= (x*y - β*z)dt + αdW_t \\\\
\\end{align}
```

with ``σ=10``, ``ρ=28``, ``β=8/3``, ``α=3.0`` and inital condition ``u₀=[1;1;1]``.
"""
prob_sde_lorenz = SDEProblem(f,σ,ones(3))

function oval2ModelExample(;largeFluctuations=false,useBigs=false,noiseLevel=1)
  #Parameters
  J1_200=3.
  J1_34=0.15
  J2_200=0.2
  J2_34=0.35
  J_2z=0.9
  J_O=0.918
  J_SO=0.5
  # J_ZEB=0.06
  J_ecad1=0.1
  J_ecad2=0.3
  J_ncad1=0.4
  J_ncad2=0.4
  J_ncad3 = 2
  J_snail0=0.6
  J_snail1=1.8
  J_zeb=3.0
  K1=1.0
  K2=1.0
  K3=1.0
  K4=1.0
  K5=1.0
  KTGF=20.
  Ks=100.
  # TGF0=0
  TGF_flg=0.
  Timescale=1000.
  dk_ZR1=0.5
  dk_ZR2=0.5
  dk_ZR3=0.5
  dk_ZR4=0.5
  dk_ZR5=0.5
  k0O=0.35
  k0_200=0.0002
  k0_34=0.001
  k0_snail=0.0005
  k0_zeb=0.003
  kO=1.2
  kOp=10.
  k_200=0.02
  k_34=0.019
  k_OT=1.1
  k_SNAIL=16.
  k_TGF=1.5
  k_ZEB=16.
  k_ecad0=5.
  k_ecad1=15.
  k_ecad2=5.
  k_ncad0=5.
  k_ncad1=2.
  k_ncad2=5.
  k_snail=0.05
  k_tgf=0.05
  k_zeb=0.06
  kdO=1.
  kd_200=0.035
  kd_34=0.035
  kd_SNAIL=1.6
  kd_SR1=0.9
  kd_TGF=0.9
  kd_ZEB=1.66
  kd_ecad=0.05
  kd_ncad=0.05
  kd_snail=0.09
  kd_tgf=0.1
  kd_tgfR=1.0
  kd_zeb=0.1
  kd_Op = 10.
  lamda1=0.5
  lamda2=0.5
  lamda3=0.5
  lamda4=0.5
  lamda5=0.5
  lamdas=0.5
  lamdatgfR=0.8
  nO=6.
  nSO=2.
  nzo=2.
  GE = 1.
  function f(t,y,dy)
    # y(1) = snailt
    # y(2) = SNAIL
    # y(3) = miR34t
    # y(4) = SR1 # abundance of SNAIL/miR34 complex
    # y(5) = zebt
    # y(6) = ZEB
    # y(7) = miR200t
    # y(8) = ZR1 # abundance of ZEB/miR200 complex with i copies of miR200 bound on the sequence of ZEB1
    # y(9) = ZR2
    # y(10) = ZR3
    # y(11) = ZR4
    # y(12) = ZR5
    # y(13) = tgft
    # y(14) = TGF
    # y(15) = tgfR # abundance of TGF/miR200 complex
    # y(16) = Ecad
    # y(17) = Ncad
    # y(18) = Ovol2
    TGF0=.5(t>100)
    #ODEs
    dy[1]=k0_snail+k_snail*(((y[14]+TGF0)/J_snail0))^2/(1+(((y[14]+TGF0)/J_snail0))^2+(y[19]/J_SO)^nSO)/(1+y[2]/J_snail1)-kd_snail*(y[1]-y[4])-kd_SR1*y[4]
    dy[2]=k_SNAIL*(y[1]-y[4])-kd_SNAIL*y[2]
    dy[3]=k0_34+k_34/(1+((y[2]/J1_34))^2+((y[6]/J2_34))^2)-kd_34*(y[3]-y[4])-kd_SR1*y[4]+lamdas*kd_SR1*y[4]
    dy[4]=Timescale*(Ks*(y[1]-y[4])*(y[3]-y[4])-y[4])
    dy[5]=k0_zeb+k_zeb*((y[2]/J_zeb))^2/(1+((y[2]/J_zeb))^2+((y[19]/J_2z))^nO)-kd_zeb*(y[5]-(5*y[8]+10*y[9]+10*y[10]+5*y[11]+y[12]))-dk_ZR1*5*y[8]-dk_ZR2*10*y[9]-dk_ZR3*10*y[10]-dk_ZR4*5*y[11]-dk_ZR5*y[12]
    dy[6]=k_ZEB*(y[5]-(5*y[8]+10*y[9]+10*y[10]+5*y[11]+y[12]))-kd_ZEB*y[6]
    dy[7]=k0_200+k_200/(1+((y[2]/J1_200))^3+((y[6]/J2_200))^2)-kd_200*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])-dk_ZR1*5*y[8]-dk_ZR2*2*10*y[9]-dk_ZR3*3*10*y[10]-dk_ZR4*4*5*y[11]-dk_ZR5*5*y[12]+lamda1*dk_ZR1*5*y[8]+lamda2*dk_ZR2*2*10*y[9]+lamda3*dk_ZR3*3*10*y[10]+lamda4*dk_ZR4*4*5*y[11]+lamda5*dk_ZR5*5*y[12]-kd_tgfR*y[15]+lamdatgfR*kd_tgfR*y[15]
    dy[8]=Timescale*(K1*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*(y[5]-(5*y[8]+10*y[9]+10*y[10]+5*y[11]+y[12]))-y[8])
    dy[9]=Timescale*(K2*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*y[8]-y[9])
    dy[10]=Timescale*(K3*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*y[9]-y[10])
    dy[11]=Timescale*(K4*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*y[10]-y[11])
    dy[12]=Timescale*(K5*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*y[11]-y[12])
    dy[13]=k_tgf-kd_tgf*(y[13]-y[15])-kd_tgfR*y[15]
    dy[14]=k_OT+k_TGF*(y[13]-y[15])-kd_TGF*y[14]
    dy[15]=Timescale*(TGF_flg+KTGF*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*(y[13]-y[15])-y[15])
    dy[16]=GE*(k_ecad0+k_ecad1/(((y[2]/J_ecad1))^2+1)+k_ecad2/(((y[6]/J_ecad2))^2+1)-kd_ecad*y[16])
    dy[17]=k_ncad0+k_ncad1*(((y[2]/J_ncad1))^2)/(((y[2]/J_ncad1))^2+1)+k_ncad2*(((y[6]/J_ncad2))^2)/(((y[6]/J_ncad2)^2+1)*(1+y[19]/J_ncad3))-kd_ncad*y[17]
    dy[18]=k0O+kO/(1+((y[6]/J_O))^nzo)-kdO*y[18]
    dy[19]=kOp*y[18]-kd_Op*y[19]
    return(dy)
  end

  function σ1(t,y,dσ)
    dσ[1] = noiseLevel*1.5y[1]
    dσ[18]= noiseLevel*6y[18]
  end

  function σ2(t,y,dσ)
    dσ[1] = 0.02y[1]
    dσ[16]= 0.02y[16]
    dσ[18]= 0.2y[18]
    dσ[17]= 0.02y[17]
  end

  if largeFluctuations
    σ = σ1
  else
    σ = σ2
  end

  u₀ = [0.128483;1.256853;0.0030203;0.0027977;0.0101511;0.0422942;0.2391346;0.0008014;0.0001464;2.67e-05;4.8e-6;9e-7;0.0619917;1.2444292;0.0486676;199.9383546;137.4267984;1.5180203;1.5180203] #Fig 9B
  if useBigs
    u₀ = big(u₀)
  end
  #u₀ =  [0.1701;1.6758;0.0027;0.0025;0.0141;0.0811;0.1642;0.0009;0.0001;0.0000;0.0000;0.0000;0.0697;1.2586;0.0478;194.2496;140.0758;1.5407;1.5407] #Fig 9A
  return(SDEProblem(f,σ,u₀))
end
