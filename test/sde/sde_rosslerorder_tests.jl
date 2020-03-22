using StochasticDiffEq, Test

SRIW1_tab = constructSRIW1()
@test minimum(checkSRIOrder(SRIW1_tab))

SRA1_tab = constructSRA1()
@test minimum(checkSRAOrder(SRA1_tab))

DRI1_tab = constructDRI1()
@test minimum(checkRIOrder(DRI1_tab))

RI1_tab = constructRI1()
@test minimum(checkRIOrder(RI1_tab))
