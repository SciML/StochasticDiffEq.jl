using StochasticDiffEq

SRIW1 = constructSRIW1()
@test minimum(checkSRIOrder(SRIW1))

SRA1 = constructSRA1()
@test minimum(checkSRAOrder(SRA1))
