"""
This function approximates solutions to the characteristic equation 
of the linearization of the DDE model for different parameter settings.
It returns the realpart of the solution with maximal realpart.

The second input params should contain parameters of the model in the 
following order:
tracing coverage, rel. freq. of test. in U2, rel. freq. of test. in Q,
strict. of isol., strict. of quarantine, tracing delay, scal. fac.
for presym. infectiousness, dur. of lat. phase, dur. of presym. phase, 
total dur. of inf. phase, basic repr. nbr.

PRCC is a boolean determining whether the code is run 
for the PRCC analysis in which case complQ is determined by the value of 
complI.
"""

import numpy as np
from cheb import cheb

def max_realpart(phi,params,PRCC):    
    trcov = params[0]
    sigma_U = params[1]
    sigma_Q = params[2]
    p_I = params[3]
    if PRCC == 0:
        p_Q = params[4]
    else: 
        if p_I < 0.5:
            p_Q = 2*p_I
        else:
            p_Q = 1
    if p_Q==3:
        print("UPS")
    tau = params[5]
    theta = params[6]
    alpha = 1/params[7]
    # baselines for dur. of presym. phase and total dur. of inf. phase are 2
    # and 9, respectively. Currently, the code is written such that only one 
    # of the two params can differ from this baseline setting b/c in our 
    # analysis when varying the total dur. of the inf. period we adjust
    # the dur. of the presym period to maintain a constant ratio between dur. 
    # of presym period to dur. total inf. period (2/9).
    if params[8] != 2 and params[9] != 9: 
        raise Warning("The code phistar.py does not allow BOTH the dur. of the"\
                      + " presym. per. AND the total dur. of the inf. per. to"\
                      + " differ from their baseline values 2 and 9")
    if params[8] == 2:
        gamma1 = 1/(2/9*params[9]) # pre-sympt. to sympt. rate
        gamma2 = 1/(7/9*params[9]) # sympt. to recov. rate
    else: 
        gamma1 = 1/params[8] # pre-sympt. to sympt. rate
        gamma2 = 1/(9-params[8]) # sympt. to recov. rate
    R0 = params[10]
    
    Ntot = 83000000 # total population size
    
    ####################################################
    ### resulting testing rates evaluated at the DFE ###
    ####################################################
    sigma_plus = 200000 # testing capacity 
    # the following three lines ensure that the decay rate of tests "sigma_minus"
    # is set such that when the system is close to the DFE approx. 85000 tests
    # are conducted per day
    ratio_pm = Ntot/85000 
    sigma_minus_fac = np.round((ratio_pm*sigma_plus-Ntot)/Ntot,3)
    sigma_minus = sigma_minus_fac*Ntot
    
    etaU1 = (sigma_plus)/(Ntot+sigma_minus) # testing rate in U1
    etaU2 = (sigma_plus*sigma_U)/(Ntot+sigma_minus) # testing rate in U2
    etaQ = (sigma_plus*sigma_Q)/(Ntot+sigma_minus) # testing rate in Q compartments
    
    ######################################################
    ### transmission rates resulting from theta and R0 ###
    ######################################################
    # this setting ensures: R0 as above and that beta1 = beta2*theta
    beta2 = R0/(theta * 1/gamma1 + 1/gamma2) # transm. rate in sym
    beta1 = beta2 * theta # transm. rate in presym
    
    beta1_t = phi*beta1 # transmission rates after accouting for level of effective contacts
    beta2_t = phi*beta2
    
    ctilde = 0.8 # baseline reported close contact rate
    ptilde2 = trcov*beta2/ctilde # infection probability for contacts made during sym
    ptilde1 = theta*ptilde2 # infection probability for contacts made during presym
    if ptilde1 > 1:
        raise Warning("The parameter constellation input to phistar.py resulted"\
                      + " in ptilde1>1!")
    if ptilde2 > 1:
        raise Warning("The parameter constellation input to phistar.py resulted"\
                      + " in ptilde2>1!")
        
    ctilde_t = phi*ctilde # reported close contact rate after accounting for 
                          # level of effective contacts
   
    # lengths of tracing intervals
    Jpre = 1/(gamma1 + etaU1) 
    Jsym = 1/(gamma2 + etaU2)
    
    # approx. time lengths different contacts are already inf. before traced
    rU1 = tau + (1/2) * Jpre # contacts from U1-index cases
    rU2_pre = tau + (1/2) * Jpre + Jsym # contacts from U2-index cases made 
                                        # during presym
    rU2_sym = tau + (1/2) * Jsym # contacts from U2-index cases made during sym
    
    # resulting tracing terms 
    c_inf_U1 = Jpre * ptilde1 * ctilde_t * etaU1    # infected contacts from U1
    c_inf_U2_U1 = Jpre * ptilde1 * ctilde_t * etaU2 # infected contacts from 
                                                  # U2-index cases made during 
                                                  # presym
    c_inf_U2_U2 = Jsym * ptilde2 * ctilde_t * etaU2 # infected contacts from 
                                                  # U2-index cases made during 
                                                  # sym
    
    # we distribute these different types of reported infected contacts on the
    # different compartments following the first order kinetics of the model 
    # and the time lengths those different contacts are already inf. before 
    # traced
    
    # contacts ending in E-comp.:
    chi_U1_E = c_inf_U1 * np.exp(-alpha*rU1) # contacts infected by U1-index cases
    chi_U2_E_pre = c_inf_U2_U1 * np.exp(-alpha*rU2_pre) # contacts infected by 
                                                     # U2-index cases during
                                                     # presym
    chi_U2_E_sym = c_inf_U2_U2 * np.exp(-alpha*rU2_sym) # contacts infected by 
                                                     # U2-index cases during 
                                                     # sym
    # contacts ending in U1-comp:
    chi_U1_U1 = c_inf_U1 * ((alpha)/(gamma1+etaU1-alpha)) * (np.exp(-alpha*rU1) \
                 - np.exp(-(gamma1+etaU1)*rU1)) # contacts infected by U1-index 
                                                # cases
    chi_U2_U1_pre  = c_inf_U2_U1 * ((alpha)/(gamma1+etaU1-alpha)) \
                      * (np.exp(-alpha*rU2_pre) \
                      - np.exp(-(gamma1+etaU1)*rU2_pre)) # contacts infected by 
                                                         # U2-index cases during
                                                         # presym
    chi_U2_U1_sym = c_inf_U2_U2 * ((alpha)/(gamma1+etaU1-alpha)) \
                    * (np.exp(-alpha*rU2_sym) \
                    - np.exp(-(gamma1+etaU1)*rU2_sym)) # contacts infected by 
                                                       # U2-index cases during 
                                                       # sym
    # contacts ending in U2-comp.:
    chi_U1_U2 = c_inf_U1 * (gamma1*alpha/(gamma1+etaU1-alpha)) \
                * ((1/(gamma2+etaU2-alpha)) * (np.exp(-alpha*rU1) \
                - np.exp(-(gamma2+etaU2)*rU1)) - (1/(gamma2+etaU2-gamma1-etaU1)) \
                * (np.exp(-(gamma1+etaU1)*rU1) \
                - np.exp(-(gamma2+etaU2)*rU1))) # contacts infected by U1-index 
                                                # cases
    chi_U2_U2_pre = c_inf_U2_U1 * (gamma1*alpha/(gamma1+etaU1-alpha)) \
                    * ((1/(gamma2+etaU2-alpha)) * (np.exp(-alpha*rU2_pre) \
                    - np.exp(-(gamma2+etaU2)*rU2_pre)) \
                    - (1/(gamma2+etaU2-gamma1-etaU1)) \
                    * (np.exp(-(gamma1+etaU1)*rU2_pre) \
                    - np.exp(-(gamma2+etaU2)*rU2_pre))) # contacts infected by 
                                                        # U2-index cases during
                                                        # presym
    chi_U2_U2_sym = c_inf_U2_U2 * (gamma1*alpha/(gamma1+etaU1-alpha)) \
                    * ((1/(gamma2+etaU2-alpha)) * (np.exp(-alpha*rU2_sym) \
                    - np.exp(-(gamma2+etaU2)*rU2_sym)) \
                    - (1/(gamma2+etaU2-gamma1-etaU1)) \
                    * (np.exp(-(gamma1+etaU1)*rU2_sym) \
                    - np.exp(-(gamma2+etaU2)*rU2_sym)))# contacts infected by 
                                                       # U2-index cases during 
                                                       # sym
    # Sum up the infected contacts of U2-index cases ending in same compartments
    chi_U2_E = chi_U2_E_pre + chi_U2_E_sym # infected contacts from U2-index 
                                           # that end up in E
    chi_U2_U1 = chi_U2_U1_pre + chi_U2_U1_sym # infected contacts from U2-index 
                                              # that end up in U1
    chi_U2_U2 = chi_U2_U2_pre + chi_U2_U2_sym # infected contacts from U2-index 
                                              # that end up in U2
    # matrices appearing in the linearization of the DDE system
    A0 = [[-alpha, 0, beta1_t, beta1_t*p_Q, beta1_t*p_I, beta2_t, beta2_t*p_Q,\
           beta2_t*p_I], 
          [0, -alpha, 0, 0, 0, 0, 0, 0],
          [alpha, 0, -etaU1-gamma1, 0, 0, 0, 0, 0], 
          [0, alpha, 0, -etaQ-gamma1, 0, 0, 0, 0], 
          [0, 0, etaU1, etaQ, -gamma1, 0, 0, 0],
          [0, 0, gamma1, 0, 0,-etaU2-gamma2, 0, 0],
          [0, 0, 0, gamma1, 0, 0, -etaQ-gamma2, 0],
          [0, 0, 0, 0, gamma1, etaU2, etaQ, -gamma2]];
        
    A1 = [[0, 0, -chi_U1_E, 0, 0, -chi_U2_E, 0, 0],
          [0, 0, chi_U1_E, 0, 0, chi_U2_E, 0, 0],
          [0, 0, -chi_U1_U1, 0, 0, -chi_U2_U1, 0, 0],
          [0, 0, chi_U1_U1, 0, 0, chi_U2_U1, 0, 0],
          [0, 0, 0, 0, 0, 0, 0, 0],
          [0, 0, -chi_U1_U2, 0, 0, -chi_U2_U2, 0, 0],
          [0, 0, chi_U1_U2, 0, 0, chi_U2_U2, 0, 0],
          [0, 0, 0, 0, 0, 0, 0, 0]];
    
    # the following lines are based on the code described in: 
    # E.Jarlebring, Some numerical methods to compute the eigenvalues of a 
    # time-delay system using Matlab, The delay e-letter, 2 (2008), 155.
    N = 10
    n = len(A0);#% Discretization nodes N and size of DDE n
    D = -cheb(N-1)*2/tau
    y = np.linalg.eigvals(np.append(np.kron(D[:-1],np.identity(n)),\
        np.append(A1,np.append(np.zeros((n,(N-2)*n)),A0,axis=1),axis=1),\
        axis=0))
    return np.max(np.real(y))
