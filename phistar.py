"""
This function calculates the level phistar of effective contacts in [0,1] at 
which the DFE of the DDE model switches from unstable to stable.

First input params should contain the parameters of the model in the 
following order:
tracing coverage, rel. freq. of test. in U2, rel. freq. of test. in Q,
strict. of isol., strict. of quarantine, tracing delay, scal. fac.
for presym. infectiousness, dur. of lat. phase, dur. of presym. phase, 
total dur. of inf. phase, basic repr. nbr.

PRCC is a boolean determining whether the code is run 
for the PRCC analysis in which case complQ is determined by the value of 
complI.
"""

from scipy.optimize import root_scalar
from max_realpart import max_realpart

def phistar(params, PRCC):
    if max_realpart(1,params,PRCC) < 0:
        return 1
    phistar = root_scalar(max_realpart, args=(params,PRCC), bracket=[0,1],\
                          method='brentq').root
    return phistar
