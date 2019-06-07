function SROCK1ConstantCache{T}(zprev) where T
  ms = SVector{10, Int}(3,5,7,10,25,50,75,100,150,200)
  mη = SVector{10, T}(2.2,12.0,13.0,14.3,20.3,27.2,32.1,36.0,42.1,46.7)
  SROCK1ConstantCache{typeof(zprev),T}(ms, mη, zprev, 1, one(T))
end
