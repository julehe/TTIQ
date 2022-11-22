"""
This function simulates the DDE model.
Vector IC: should contain the initial conditions for:
S, E, QE, U1, QU1, I1, U2, QU2, I2, R, Icum 
The last compartment is a dumy to count cumulative nbr. of positively 
tested individuals.

Function history: history function to initialize the DDE

Scalar t_interv: possible change point for parameter values

Vector para_cap: should contain values for TTIQ capacity paramaters in the
folowing order:
sigma_plus, p, Omega

Vector para: should contain parameters of the model in the 
following order (values for 0<=t<t_interv):
phi, tracing coverage, rel. freq. of test. in U2, rel. freq. of test. in Q,
strict. of isolation, strict. of quarantine, tracing delay.

para_hist: same as para but values for t<0.

para_interv: same as para but for t>=t_interv

Boolean Fig11: should only be set to 1 if tracing delay is an integer. Gives 
additional outputs (explained below).
"""

import numpy as np
from ddeint import ddeint

def model(steps, duration, IC, history, para_cap, para, para_hist, para_interv,\
          t_interv, Fig11):
    T = np.linspace(0,duration,steps*duration+1)
    Ntot = np.sum(IC) 
    # baseline parameter values
    alpha = 1/3.5 # latent to infectious
    gamma1 = 1/2 # presym to sym
    gamma2 = 1/7 # sym to recovered
    R0 = 3.3 # basic reproduction number
    theta = 1.5 # scaling factor for presym. transmission rate
    tau = para[6] #tracing delay
    ######################################################
    ### transmission rates resulting from theta and R0 ###
    ######################################################
    beta2 = R0/(theta * 1/gamma1 + 1/gamma2) # transm. rate in sym
    beta1 = beta2 * theta # transm. rate in presym.
    def phi(t): # time-dep. level of effective contacts
        if t<0:
            return para_hist[0]
        elif t<t_interv:
            return para[0]
        else:
            return para_interv[0]
    # transmission rates after accounting for phi
    def beta1_t(t):
        return phi(t)*beta1
    def beta2_t(t):
        return phi(t)*beta2
    ###############
    ### Tracing ###
    ###############
    Omega = para_cap[0] # capacity
    p = para_cap[1] # efficiency constant
    def trcov_t(t): # time dep. tracing-coverage
        if t<0:
            return para_hist[1]
        elif t<t_interv:
            return para[1]
        else:
            return para_interv[1]
    ctilde = 0.8 # baseline reported close contact rate
    ###############
    ### Testing ###
    ###############
    sigma_plus = para_cap[2] # testing capacity 
    # the following three lines ensure that the decay rate of tests "sigma_minus"
    # is set such that when the system is close to the DFE approx. 85000 tests
    # are conducted per day
    ratio_pm = Ntot/85000 
    sigma_minus_fac = np.round((ratio_pm*sigma_plus-Ntot)/Ntot,3)
    sigma_minus = sigma_minus_fac*Ntot
    def sigma_U_t(t): # time-dep. rel. freq. of testing in U2
        if t<0:
            return para_hist[2]
        elif t<t_interv:
            return para[2]
        else:
            return para_interv[2]
    def sigma_Q_t(t): # time-dep. rel. freq. of testing in Q
        if t<0:
            return para_hist[3]
        elif t<t_interv:
            return para[3]
        else:
            return para_interv[3]
    def p_Q_t(t): # time-dep. strictness of quarantine
        if t<0:
            return para_hist[5]
        elif t<t_interv:
            return para[5]
        else:
            return para_interv[5]
    def p_I_t(t): # time-dep. strictness of isolation
        if t<0:
            return para_hist[4]
        elif t<t_interv:
            return para[4]
        else:
            return para_interv[4]
    # differential equations
    def rhs(X, t, tau, caps):
        S, E, QE, U1, QU1, I1, U2, QU2, I2, R, Icum = X(t)
        St, Et, QEt, U1t, QU1t, I1t, U2t, QU2t, I2t, Rt, Icumt = X(t-tau)
        # current/past values of the time-dep. paras
        betaU1 = beta1_t(t)
        betaU2 = beta2_t(t)
        sigma_U = sigma_U_t(t)
        sigma_Ut = sigma_U_t(t-tau)
        sigma_Q = sigma_Q_t(t)
        sigma_Qt = sigma_Q_t(t-tau)
        trcov = trcov_t(t) 
        p_Q = p_Q_t(t)
        p_I = p_I_t(t)
        # calculate current and past testing terms accounting for system's state
        Nenner = sigma_minus + S + E + sigma_Q*QE + U1 + sigma_Q*QU1 + I1\
                            + sigma_U*U2 + sigma_Q*QU2 + I2 + R
        Nennert = sigma_minus + St + Et + sigma_Qt*QEt + U1t + sigma_Qt*QU1t + I1t\
                            + sigma_Ut*U2t + sigma_Qt*QU2t + I2t + Rt
        # correction factor for testing rates (accounting for limited capacity)           
        cor = sigma_plus/Nenner
        cort = sigma_plus/Nennert
        # time lagged incidence of new cases from U
        Testt1 = cort*U1t
        Testt2 = cort*sigma_Ut*U2t
        # lengths tracing intervals
        J = 9
        Jpre = 1/(gamma1 + cort)
        Jsym = 1/(gamma2 + cort*sigma_Ut)
        # approx. time lengths different contacts are already inf. before traced
        rU1 = tau + (1/2) * Jpre # contacts from U1-index cases
        rU2_pre = tau + (1/2) * Jpre + Jsym # contacts from U2-index cases made 
                                            # during presym
        rU2_sym = tau + (1/2) * Jsym # contacts from U2-index cases made during sym
        
        ctildet = phi(t-tau)*ctilde # reported close contact rate after accounting
                                    # for level of effective contacts
        # infection probabilities for close contacts (as a result of the tracing coverage)
        ptilde2 = trcov*beta2/ctilde # infection probability for contacts made during sym
        ptilde1 = theta*ptilde2 # infection probability for contacts made during presym
        if ptilde1 > 1:
            raise Warning("The parameter input to model.py resulted in ptilde1>1!")
        if ptilde2 > 1:
            raise Warning("The parameter input to model.py resulted in ptilde2>1!")
        # resulting tracing terms              
        c_pot_U1 = J * ctildet * Testt1
        c_pot_U2 = J * ctildet * Testt2
        c_pot = c_pot_U1 + c_pot_U2
        c_act_U1 = c_pot_U1 * Omega/(Omega**p + c_pot**p)**(1/p)
        c_act_U2 = c_pot_U2 * Omega/(Omega**p + c_pot**p)**(1/p)
        c_act_inf_U1 = Jpre/J * ptilde1 * (St/Ntot) * c_act_U1
        c_act_inf_U2_U1 = Jpre/J * ptilde1 * (St/Ntot) * c_act_U2
        c_act_inf_U2_U2 = Jsym/J * ptilde2 * (St/Ntot) * c_act_U2
        # we distribute these different types of reported infected contacts on the
        # different compartments following the first order kinetics of the model 
        # and the time lengths those different contacts are already inf. before 
        # traced
        # contacts ending in E-comp.:
        Tr_U1_E = c_act_inf_U1 * np.exp(-alpha*rU1) # contacts infected by U1-index cases
        Tr_U2_E_pre = c_act_inf_U2_U1 * np.exp(-alpha*rU2_pre) # contacts infected by 
                                                         # U2-index cases during
                                                         # presym
        Tr_U2_E_sym = c_act_inf_U2_U2 * np.exp(-alpha*rU2_sym) # contacts infected by 
                                                         # U2-index cases during 
                                                         # sym
        
        # contacts ending in U1-comp:
        Tr_U1_U1 = c_act_inf_U1 * ((alpha)/(gamma1+cort-alpha)) \
                   * (np.exp(-alpha*rU1) - np.exp(-(gamma1+cort)*rU1))
                   # contacts infected by U1-index cases
        Tr_U2_U1_pre  = c_act_inf_U2_U1 \
                        * ((alpha)/(gamma1+cort-alpha)) \
                        * (np.exp(-alpha*rU2_pre) \
                        - np.exp(-(gamma1+cort)*rU2_pre))
                        # contacts infected by U2-index cases during presym
        Tr_U2_U1_sym = c_act_inf_U2_U2 * ((alpha)/(gamma1+cort-alpha))\
                       * (np.exp(-alpha*rU2_sym) \
                       - np.exp(-(gamma1+cort)*rU2_sym))
                       # contacts infected by U2-index cases during sym
        # contacts ending in U2-comp.:
        Tr_U1_U2 = c_act_inf_U1 * (gamma1*alpha/(gamma1+cort-alpha)) \
                   * ((1/(gamma2+cort*sigma_Ut-alpha)) * (np.exp(-alpha*rU1) \
                   - np.exp(-(gamma2+cort*sigma_Ut)*rU1)) \
                   - (1/(gamma2+cort*sigma_Ut-gamma1-cort))\
                   * (np.exp(-(gamma1+cort)*rU1)\
                   - np.exp(-(gamma2+cort*sigma_Ut)*rU1)))
                   # contacts infected by U1-index cases
        Tr_U2_U2_pre = c_act_inf_U2_U1 \
                       * (gamma1*alpha/(gamma1+cort-alpha)) \
                       * ((1/(gamma2+cort*sigma_Ut-alpha)) \
                       * (np.exp(-alpha*rU2_pre) \
                       - np.exp(-(gamma2+cort*sigma_Ut)*rU2_pre)) \
                       - (1/(gamma2+cort*sigma_Ut-gamma1-cort)) \
                       * (np.exp(-(gamma1+cort)*rU2_pre) \
                       - np.exp(-(gamma2+cort*sigma_Ut)*rU2_pre)))
                       # contacts infected by U2-index cases during presym
        Tr_U2_U2_sym = c_act_inf_U2_U2 * (gamma1*alpha/(gamma1+cort-alpha))\
                        * ((1/(gamma2+cort*sigma_Ut-alpha)) \
                        * (np.exp(-alpha*rU2_sym) \
                        - np.exp(-(gamma2+cort*sigma_Ut)*rU2_sym)) \
                        - (1/(gamma2+cort*sigma_Ut-gamma1-cort)) \
                        * (np.exp(-(gamma1+cort)*rU2_sym) \
                        - np.exp(-(gamma2+cort*sigma_Ut)*rU2_sym)))
                        # contacts infected by U2-index cases during sym
        # Sum up the infected contacts of U2-index cases ending in same compartments
        Tr_E = Tr_U1_E + Tr_U2_E_pre + Tr_U2_E_sym
        Tr_U1 = Tr_U1_U1 + Tr_U2_U1_pre + Tr_U2_U1_sym
        Tr_U2 = Tr_U1_U2 + Tr_U2_U2_pre + Tr_U2_U2_sym
        # force of infection
        lam = (betaU1*U1 + p_Q*betaU1*QU1 + p_I*betaU1*I1 + betaU2*U2 \
               + p_Q*betaU2*QU2 + p_I*betaU2*I2)/Ntot 
        # right-hand side
        return np.array([-lam*S,                                                 #S
                          lam*S - alpha*E - Tr_E,                                #E
                          Tr_E - alpha*QE,                                       #QE
                          alpha*E - (cor + gamma1)*U1 - Tr_U1,                   #U1
                          Tr_U1 + alpha*QE - gamma1*QU1 - cor*sigma_Q*QU1,       #QU1
                          cor*(U1 + sigma_Q*QU1) - gamma1*I1,                    #I1
                          gamma1*U1 - cor*sigma_U*U2 - Tr_U2 - gamma2*U2,        #U2  
                          Tr_U2 + gamma1*QU1 - cor*sigma_Q*QU2 - gamma2*QU2,     #QU2
                          cor*(sigma_U*U2 + sigma_Q*QU2) + gamma1*I1 - gamma2*I2,#I2
                          gamma2* (U2 + QU2 + I2),                               #R
                          cor*(sigma_U*U2 + sigma_Q*QU2 + U1 + sigma_Q*QU1)])    #Icum

    v = ddeint(rhs, history, T, fargs=(tau,1,))
    s = np.array(v[:,0])
    e = np.array(v[:,1])
    qe = np.array(v[:,2])
    u1 = np.array(v[:,3])
    qu1 = np.array(v[:,4])
    i1 = np.array(v[:,5])
    u2 = np.array(v[:,6])
    qu2 = np.array(v[:,7])
    i2 = np.array(v[:,8])
    r = np.array(v[:,9])
    # reconstruct times series for incidence term lam*S
    lamS = []
    for j in range(0,len(s)):
        lamS.append(s[j]*(beta1_t(T[j])*u1[j] + p_Q_t(T[j])*beta1_t(T[j])*qu1[j] \
                    + p_I_t(T[j])*beta1_t(T[j])*i1[j] + beta2_t(T[j])*u2[j] \
                    + p_Q_t(T[j])*beta2_t(T[j])*qu2[j] \
                    + p_I_t(T[j])*beta2_t(T[j])*i2[j])/Ntot)
    # In the script Fig12_model_sim, where tau=2, we plot the tracing effciency, 
    # test positive rate, and detection ratio. 
    # The code in the following if-block reconstructs the corresponding time 
    # series. Currently, this only works for integer tracing delays tau. 
    # Since other scripts consider scenarios with non-integer tau, we condition
    # this part to the boolean Fig12. One could use interpolation of the solution v
    # to make it run for non-integer tau at the cost of increased computation time.
    if Fig12 == True:
        cor = []    
        eta_U1 = []
        eta_U2 = []
        eta_QU1 = []
        eta_QU2 = [] 
        scalar = []
        for t in range(0,len(T)):
            sigma_U = sigma_U_t(T[t])
            sigma_Q = sigma_Q_t(T[t])
            Nenner = sigma_minus + s[t] + e[t] + sigma_Q*qe[t]\
                        + u1[t] + sigma_Q*qu1[t] + i1[t]\
                          +  sigma_U*u2[t] + sigma_Q*qu2[t] + i2[t] + r[t]
            scalar.append(Nenner - sigma_minus)  
            cor.append(sigma_plus/Nenner)
            eta_U1.append(cor[-1])
            eta_U2.append(cor[-1]*sigma_U)
            eta_QU1.append(cor[-1]*sigma_Q)
            eta_QU2.append(cor[-1]*sigma_Q)
        # reconstruct time lagged testing rates and incidences   
        corst = []    
        eta_U1t = []
        eta_U2t = []
        inc1t = []
        inc2t = [] 
        for t in range(0,len(T)):
            if t<tau*steps:
                sigma_Ut = sigma_U_t(T[t]-tau)
                sigma_Qt = sigma_Q_t(T[t]-tau)
                Nenner0 = sigma_minus + IC[0] + IC[1] + sigma_Qt*IC[2]\
                        + IC[3] + sigma_Qt*IC[4] + IC[5]\
                        + sigma_Ut*IC[6] + sigma_Qt*IC[7] + IC[8] + IC[9];
                corst.append(sigma_plus/Nenner0)
                eta_U1t.append(corst[-1])
                eta_U2t.append(corst[-1]*sigma_Ut)
                inc1t.append(eta_U1t[-1]*IC[3])
                inc2t.append(eta_U2t[-1]*IC[6])
            else:
                sigma_Ut = sigma_U_t(T[t]-tau)
                sigma_Qt = sigma_Q_t(T[t]-tau)
                Nenner = sigma_minus + s[t-(tau*steps)] + e[t-(tau*steps)] \
                         + sigma_Qt*qe[t-(tau*steps)] + u1[t-(tau*steps)] \
                         + sigma_Qt*qu1[t-(tau*steps)] + i1[t-(tau*steps)]\
                         + sigma_Ut*u2[t-(tau*steps)] + sigma_Qt*qu2[t-(tau*steps)]\
                         + i2[t-(tau*steps)] + r[t-(tau*steps)]     
                corst.append(sigma_plus/Nenner)
                eta_U1t.append(corst[-1])
                eta_U2t.append(corst[-1]*sigma_Ut)
                inc1t.append(eta_U1t[-1]*u1[t-(tau*steps)])
                inc2t.append(eta_U2t[-1]*u2[t-(tau*steps)])
        inc1t = np.array(inc1t)
        inc2t = np.array(inc2t)
        scalar = np.array(scalar)
        treff = []
        for t in range(0,len(T)):
            treff.append(Omega/(Omega**p+(9*ctilde*phi(T[t]-tau)*(inc1t[t]\
                                                    +inc2t[t]))**p)**(1/p))
        tests_used = (scalar/(np.ones(len(scalar))*sigma_minus+scalar))*sigma_plus
        test_pos = (eta_U1*u1 + eta_U2*u2 + eta_QU1*qu1 + eta_QU2*qu2)/tests_used
        det_ratio = eta_U1/(eta_U1+np.ones(len(eta_U1))*gamma1) \
                    + (np.ones(len(eta_U1))*gamma1/(eta_U1+np.ones(len(eta_U1))\
                    *gamma1)) * eta_U2/(eta_U2+np.ones(len(eta_U2))*gamma2)
    else:
        treff = []
        test_pos = []
        det_ratio = []
    return (v,treff,test_pos,det_ratio,lamS)




