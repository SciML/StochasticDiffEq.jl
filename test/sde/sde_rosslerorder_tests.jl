using StochasticDiffEq

SRIW1_tab = constructSRIW1()
@test minimum(checkSRIOrder(SRIW1_tab))

SRA1_tab = constructSRA1()
@test minimum(checkSRAOrder(SRA1_tab))
