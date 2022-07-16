# Figure 13: phibar along simulation of DDE model in baseline setting but
# for different TTIQ capacity paramaters.

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import string
from model import model

#costum colors
My_Amber = "#E69F00"
My_Blue = "#56B4E9"
My_Green = "#009E73"
My_Violet = '#661D98'

color_list = [My_Amber, My_Blue, My_Green, My_Violet]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)
# baseline parameter setting 
# level of effective contacts
phi_b = 0.6
# tracing coverage
trcov_b = 0.65
# rel. freq. of testing in U
sigma_U_b = 93
# rel. freq. of testing in Q
sigma_Q_b = 300
# strictness of isolation 
p_I_b = 0.1
# strictness of quarantine
p_Q_b = 0.2
# tracing delay
tau_b = 2
para = [phi_b, trcov_b, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b,tau_b]
para_hist = para
para_interv = para
# running simulations first on total time horizon with baseline phi
# saves computing time since later simulations can use interpolated solution
# as history and do not have to run whole simulation again.
steps = 32
duration = 250
T = np.linspace(0,duration,steps*duration+1);
IC = [82975287, 2564, 131, 1301, 86, 52, 2173, 207, 1666, 16533, 0]
def history_base(t):
    return IC
# compare different testing capacities
sol_sigma_plus_low, *_ \
= model(steps, duration, IC, history_base,[40000,2,150000], para, para_hist,\
        para_interv, 0, 0)
sol_sigma_plus_high, *_ \
= model(steps, duration, IC, history_base,[40000,2,250000], para, para_hist,\
        para_interv, 0, 0)

s_sigma_plus_low = interp1d(T,np.array(sol_sigma_plus_low[:,0]),\
                            kind='linear')
e_sigma_plus_low = interp1d(T,np.array(sol_sigma_plus_low[:,1]),\
                            kind='linear')
qe_sigma_plus_low = interp1d(T,np.array(sol_sigma_plus_low[:,2]),\
                             kind='linear')
u1_sigma_plus_low = interp1d(T,np.array(sol_sigma_plus_low[:,3]),\
                             kind='linear')
qu1_sigma_plus_low = interp1d(T,np.array(sol_sigma_plus_low[:,4]),\
                              kind='linear')
i1_sigma_plus_low = interp1d(T,np.array(sol_sigma_plus_low[:,5]),\
                             kind='linear')
u2_sigma_plus_low = interp1d(T,np.array(sol_sigma_plus_low[:,6]),\
                             kind='linear')
qu2_sigma_plus_low = interp1d(T,np.array(sol_sigma_plus_low[:,7]),\
                              kind='linear')
i2_sigma_plus_low = interp1d(T,np.array(sol_sigma_plus_low[:,8]),\
                             kind='linear')
r_sigma_plus_low = interp1d(T,np.array(sol_sigma_plus_low[:,9]),\
                            kind='linear')

s_sigma_plus_high = interp1d(T,np.array(sol_sigma_plus_high[:,0]),\
                            kind='linear')
e_sigma_plus_high = interp1d(T,np.array(sol_sigma_plus_high[:,1]),\
                            kind='linear')
qe_sigma_plus_high = interp1d(T,np.array(sol_sigma_plus_high[:,2]),\
                            kind='linear')
u1_sigma_plus_high = interp1d(T,np.array(sol_sigma_plus_high[:,3]),\
                            kind='linear')
qu1_sigma_plus_high = interp1d(T,np.array(sol_sigma_plus_high[:,4]),\
                            kind='linear')
i1_sigma_plus_high = interp1d(T,np.array(sol_sigma_plus_high[:,5]),\
                            kind='linear')
u2_sigma_plus_high = interp1d(T,np.array(sol_sigma_plus_high[:,6]),\
                            kind='linear')
qu2_sigma_plus_high = interp1d(T,np.array(sol_sigma_plus_high[:,7]),\
                            kind='linear')
i2_sigma_plus_high = interp1d(T,np.array(sol_sigma_plus_high[:,8]),\
                            kind='linear')
r_sigma_plus_high = interp1d(T,np.array(sol_sigma_plus_high[:,9]),\
                            kind='linear')

# compare different tracing capacities
sol_Omega_low_p_medium, *_\
= model(steps, duration, IC, history_base,[20000,2,200000], para, para_hist,\
        para_interv, 0, 0)
sol_Omega_medium_p_medium, *_\
= model(steps, duration, IC, history_base,[40000,2,200000], para, para_hist,\
        para_interv, 0, 0)
sol_Omega_high_p_medium, *_\
= model(steps, duration, IC, history_base,[60000,2,200000], para, para_hist,\
        para_interv, 0, 0)

s_Omega_low_p_medium = interp1d(T,np.array(sol_Omega_low_p_medium[:,0]),\
                            kind='linear')
e_Omega_low_p_medium = interp1d(T,np.array(sol_Omega_low_p_medium[:,1]),\
                            kind='linear')
qe_Omega_low_p_medium = interp1d(T,np.array(sol_Omega_low_p_medium[:,2]),\
                            kind='linear')
u1_Omega_low_p_medium = interp1d(T,np.array(sol_Omega_low_p_medium[:,3]),\
                            kind='linear')
qu1_Omega_low_p_medium = interp1d(T,np.array(sol_Omega_low_p_medium[:,4]),\
                            kind='linear')
i1_Omega_low_p_medium = interp1d(T,np.array(sol_Omega_low_p_medium[:,5]),\
                            kind='linear')
u2_Omega_low_p_medium = interp1d(T,np.array(sol_Omega_low_p_medium[:,6]),\
                            kind='linear')
qu2_Omega_low_p_medium = interp1d(T,np.array(sol_Omega_low_p_medium[:,7]),\
                            kind='linear')
i2_Omega_low_p_medium = interp1d(T,np.array(sol_Omega_low_p_medium[:,8]),\
                            kind='linear')
r_Omega_low_p_medium = interp1d(T,np.array(sol_Omega_low_p_medium[:,9]),\
                            kind='linear')

s_Omega_medium_p_medium = interp1d(T,np.array(sol_Omega_medium_p_medium[:,0]),\
                            kind='linear')
e_Omega_medium_p_medium = interp1d(T,np.array(sol_Omega_medium_p_medium[:,1]),\
                            kind='linear')
qe_Omega_medium_p_medium = interp1d(T,np.array(sol_Omega_medium_p_medium[:,2]),\
                            kind='linear')
u1_Omega_medium_p_medium = interp1d(T,np.array(sol_Omega_medium_p_medium[:,3]),\
                            kind='linear')
qu1_Omega_medium_p_medium = interp1d(T,np.array(sol_Omega_medium_p_medium[:,4]),\
                            kind='linear')
i1_Omega_medium_p_medium = interp1d(T,np.array(sol_Omega_medium_p_medium[:,5]),\
                            kind='linear')
u2_Omega_medium_p_medium = interp1d(T,np.array(sol_Omega_medium_p_medium[:,6]),\
                            kind='linear')
qu2_Omega_medium_p_medium = interp1d(T,np.array(sol_Omega_medium_p_medium[:,7]),\
                            kind='linear')
i2_Omega_medium_p_medium = interp1d(T,np.array(sol_Omega_medium_p_medium[:,8]),\
                            kind='linear')
r_Omega_medium_p_medium = interp1d(T,np.array(sol_Omega_medium_p_medium[:,9]),\
                            kind='linear')

s_Omega_high_p_medium = interp1d(T,np.array(sol_Omega_high_p_medium[:,0]),\
                            kind='linear')
e_Omega_high_p_medium = interp1d(T,np.array(sol_Omega_high_p_medium[:,1]),\
                            kind='linear')
qe_Omega_high_p_medium = interp1d(T,np.array(sol_Omega_high_p_medium[:,2]),\
                            kind='linear')
u1_Omega_high_p_medium = interp1d(T,np.array(sol_Omega_high_p_medium[:,3]),\
                            kind='linear')
qu1_Omega_high_p_medium = interp1d(T,np.array(sol_Omega_high_p_medium[:,4]),\
                            kind='linear')
i1_Omega_high_p_medium = interp1d(T,np.array(sol_Omega_high_p_medium[:,5]),\
                            kind='linear')
u2_Omega_high_p_medium = interp1d(T,np.array(sol_Omega_high_p_medium[:,6]),\
                            kind='linear')
qu2_Omega_high_p_medium = interp1d(T,np.array(sol_Omega_high_p_medium[:,7]),\
                            kind='linear')
i2_Omega_high_p_medium = interp1d(T,np.array(sol_Omega_high_p_medium[:,8]),\
                            kind='linear')
r_Omega_high_p_medium = interp1d(T,np.array(sol_Omega_high_p_medium[:,9]),\
                            kind='linear')

# compare different p
sol_Omega_medium_p_low, *_\
= model(steps, duration, IC, history_base,[40000,1,200000], para, para_hist,\
        para_interv, 0, 0)
sol_Omega_medium_p_high, *_\
= model(steps, duration, IC, history_base,[40000,10,200000], para, para_hist,\
        para_interv, 0, 0)

s_Omega_medium_p_low = interp1d(T,np.array(sol_Omega_medium_p_low[:,0]),\
                            kind='linear')
e_Omega_medium_p_low = interp1d(T,np.array(sol_Omega_medium_p_low[:,1]),\
                            kind='linear')
qe_Omega_medium_p_low = interp1d(T,np.array(sol_Omega_medium_p_low[:,2]),\
                            kind='linear')
u1_Omega_medium_p_low = interp1d(T,np.array(sol_Omega_medium_p_low[:,3]),\
                            kind='linear')
qu1_Omega_medium_p_low = interp1d(T,np.array(sol_Omega_medium_p_low[:,4]),\
                            kind='linear')
i1_Omega_medium_p_low = interp1d(T,np.array(sol_Omega_medium_p_low[:,5]),\
                            kind='linear')
u2_Omega_medium_p_low = interp1d(T,np.array(sol_Omega_medium_p_low[:,6]),\
                            kind='linear')
qu2_Omega_medium_p_low = interp1d(T,np.array(sol_Omega_medium_p_low[:,7]),\
                            kind='linear')
i2_Omega_medium_p_low = interp1d(T,np.array(sol_Omega_medium_p_low[:,8]),\
                            kind='linear')
r_Omega_medium_p_low = interp1d(T,np.array(sol_Omega_medium_p_low[:,9]),\
                            kind='linear')

s_Omega_medium_p_high = interp1d(T,np.array(sol_Omega_medium_p_high[:,0]),\
                            kind='linear')
e_Omega_medium_p_high = interp1d(T,np.array(sol_Omega_medium_p_high[:,1]),\
                            kind='linear')
qe_Omega_medium_p_high = interp1d(T,np.array(sol_Omega_medium_p_high[:,2]),\
                            kind='linear')
u1_Omega_medium_p_high = interp1d(T,np.array(sol_Omega_medium_p_high[:,3]),\
                            kind='linear')
qu1_Omega_medium_p_high = interp1d(T,np.array(sol_Omega_medium_p_high[:,4]),\
                            kind='linear')
i1_Omega_medium_p_high = interp1d(T,np.array(sol_Omega_medium_p_high[:,5]),\
                            kind='linear')
u2_Omega_medium_p_high = interp1d(T,np.array(sol_Omega_medium_p_high[:,6]),\
                            kind='linear')
qu2_Omega_medium_p_high = interp1d(T,np.array(sol_Omega_medium_p_high[:,7]),\
                            kind='linear')
i2_Omega_medium_p_high = interp1d(T,np.array(sol_Omega_medium_p_high[:,8]),\
                            kind='linear')
r_Omega_medium_p_high = interp1d(T,np.array(sol_Omega_medium_p_high[:,9]),\
                            kind='linear')

# cap vs. p
sol_Omega_low_p_high, *_\
= model(steps, duration, IC, history_base,[20000,10,200000], para, para_hist,\
        para_interv, 0, 0)
sol_Omega_high_p_low, *_\
= model(steps, duration, IC, history_base,[60000,1,200000], para, para_hist,\
        para_interv, 0, 0)

s_Omega_low_p_high = interp1d(T,np.array(sol_Omega_low_p_high[:,0]),\
                            kind='linear')
e_Omega_low_p_high = interp1d(T,np.array(sol_Omega_low_p_high[:,1]),\
                            kind='linear')
qe_Omega_low_p_high = interp1d(T,np.array(sol_Omega_low_p_high[:,2]),\
                            kind='linear')
u1_Omega_low_p_high = interp1d(T,np.array(sol_Omega_low_p_high[:,3]),\
                            kind='linear')
qu1_Omega_low_p_high = interp1d(T,np.array(sol_Omega_low_p_high[:,4]),\
                            kind='linear')
i1_Omega_low_p_high = interp1d(T,np.array(sol_Omega_low_p_high[:,5]),\
                            kind='linear')
u2_Omega_low_p_high = interp1d(T,np.array(sol_Omega_low_p_high[:,6]),\
                            kind='linear')
qu2_Omega_low_p_high = interp1d(T,np.array(sol_Omega_low_p_high[:,7]),\
                            kind='linear')
i2_Omega_low_p_high = interp1d(T,np.array(sol_Omega_low_p_high[:,8]),\
                            kind='linear')
r_Omega_low_p_high = interp1d(T,np.array(sol_Omega_low_p_high[:,9]),\
                            kind='linear')

s_Omega_high_p_low = interp1d(T,np.array(sol_Omega_high_p_low[:,0]),\
                            kind='linear')
e_Omega_high_p_low = interp1d(T,np.array(sol_Omega_high_p_low[:,1]),\
                            kind='linear')
qe_Omega_high_p_low = interp1d(T,np.array(sol_Omega_high_p_low[:,2]),\
                            kind='linear')
u1_Omega_high_p_low = interp1d(T,np.array(sol_Omega_high_p_low[:,3]),\
                            kind='linear')
qu1_Omega_high_p_low = interp1d(T,np.array(sol_Omega_high_p_low[:,4]),\
                            kind='linear')
i1_Omega_high_p_low = interp1d(T,np.array(sol_Omega_high_p_low[:,5]),\
                            kind='linear')
u2_Omega_high_p_low = interp1d(T,np.array(sol_Omega_high_p_low[:,6]),\
                            kind='linear')
qu2_Omega_high_p_low = interp1d(T,np.array(sol_Omega_high_p_low[:,7]),\
                            kind='linear')
i2_Omega_high_p_low = interp1d(T,np.array(sol_Omega_high_p_low[:,8]),\
                            kind='linear')
r_Omega_high_p_low = interp1d(T,np.array(sol_Omega_high_p_low[:,9]),\
                            kind='linear')

# function calculating phibar for specific timing of intervention + TTIQ 
# capacity setting
def phibar(timing,Omega,p,sigma_plus):
    timing = int(timing)
    upper = 1
    lower = 0
    middle = lower + (upper - lower)/2
    while True:
        steps = 32
        duration = 56
        IC = [82975287, 2564, 131, 1301, 86, 52, 2173, 207, 1666, 16533, 0]
        phi = middle
        para = [phi, trcov_b, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b,tau_b]
        para_hist = [phi_b, trcov_b, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b,tau_b]
        para_interv = para
        def history(t):
            t=timing+t
            #print(t)
            if t<0:
                return IC
            else:
                if Omega == 20000 and p == 2:
                    return [s_Omega_low_p_medium(t),e_Omega_low_p_medium(t),\
                            qe_Omega_low_p_medium(t),u1_Omega_low_p_medium(t),\
                            qu1_Omega_low_p_medium(t),i1_Omega_low_p_medium(t),\
                            u2_Omega_low_p_medium(t),qu2_Omega_low_p_medium(t),\
                            i2_Omega_low_p_medium(t),r_Omega_low_p_medium(t),0]
                elif Omega == 40000:
                    if p==1:
                        return [s_Omega_medium_p_low(t),e_Omega_medium_p_low(t),\
                                qe_Omega_medium_p_low(t),u1_Omega_medium_p_low(t),\
                                qu1_Omega_medium_p_low(t),i1_Omega_medium_p_low(t),\
                                u2_Omega_medium_p_low(t),qu2_Omega_medium_p_low(t),\
                                i2_Omega_medium_p_low(t),r_Omega_medium_p_low(t),0]
                    elif p==2:
                        if sigma_plus == 150000:
                            return [s_sigma_plus_low(t),e_sigma_plus_low(t),\
                                    qe_sigma_plus_low(t),u1_sigma_plus_low(t),\
                                    qu1_sigma_plus_low(t),i1_sigma_plus_low(t),\
                                    u2_sigma_plus_low(t),qu2_sigma_plus_low(t),\
                                    i2_sigma_plus_low(t),r_sigma_plus_low(t),0]
                        elif sigma_plus == 200000:
                            return [s_Omega_medium_p_medium(t),\
                                    e_Omega_medium_p_medium(t),\
                                    qe_Omega_medium_p_medium(t),\
                                    u1_Omega_medium_p_medium(t),\
                                    qu1_Omega_medium_p_medium(t),\
                                    i1_Omega_medium_p_medium(t),\
                                    u2_Omega_medium_p_medium(t),\
                                    qu2_Omega_medium_p_medium(t),\
                                    i2_Omega_medium_p_medium(t),\
                                    r_Omega_medium_p_medium(t),0]
                        elif sigma_plus == 250000:
                            return [s_sigma_plus_high(t),e_sigma_plus_high(t),\
                                    qe_sigma_plus_high(t),u1_sigma_plus_high(t),\
                                    qu1_sigma_plus_high(t),i1_sigma_plus_high(t),\
                                    u2_sigma_plus_high(t),qu2_sigma_plus_high(t),\
                                    i2_sigma_plus_high(t),r_sigma_plus_high(t),0]
                    else:
                        return [s_Omega_medium_p_high(t),e_Omega_medium_p_high(t),\
                                qe_Omega_medium_p_high(t),u1_Omega_medium_p_high(t),\
                                qu1_Omega_medium_p_high(t),i1_Omega_medium_p_high(t),\
                                u2_Omega_medium_p_high(t),qu2_Omega_medium_p_high(t),\
                                i2_Omega_medium_p_high(t),r_Omega_medium_p_high(t),0]
                elif Omega == 60000 and p == 2:
                    return [s_Omega_high_p_medium(t),e_Omega_high_p_medium(t),\
                            qe_Omega_high_p_medium(t),u1_Omega_high_p_medium(t),\
                            qu1_Omega_high_p_medium(t),i1_Omega_high_p_medium(t),\
                            u2_Omega_high_p_medium(t),qu2_Omega_high_p_medium(t),\
                            i2_Omega_high_p_medium(t),r_Omega_high_p_medium(t),0]
                elif Omega == 20000 and p == 10:
                    return [s_Omega_low_p_high(t),e_Omega_low_p_high(t),\
                            qe_Omega_low_p_high(t),u1_Omega_low_p_high(t),\
                            qu1_Omega_low_p_high(t),i1_Omega_low_p_high(t),\
                            u2_Omega_low_p_high(t),qu2_Omega_low_p_high(t),\
                            i2_Omega_low_p_high(t),r_Omega_low_p_high(t),0]
                else: 
                    return [s_Omega_high_p_low(t),e_Omega_high_p_low(t),\
                            qe_Omega_high_p_low(t),u1_Omega_high_p_low(t),\
                            qu1_Omega_high_p_low(t),i1_Omega_high_p_low(t),\
                            u2_Omega_high_p_low(t),qu2_Omega_high_p_low(t),\
                            i2_Omega_high_p_low(t),r_Omega_high_p_low(t),0]
        
        *_, lamS\
        = model(steps, duration, IC, history, [Omega,p,sigma_plus], para,\
                para_hist, para_interv, 0, 0)
        if (lamS)[14*steps] >= max((lamS)[14*steps:]) or (lamS)[0] >= max((lamS)[0:]):
            if upper - lower < 0.000001:
                break
            else:
                upper = upper
                lower = middle
                middle = lower + (upper - lower)/2
                continue
        else:
            lower = lower
            upper = middle
            middle = lower + (upper - lower)/2
            continue
    return middle

# set stepsize and time interval on which to calculate phibar 
n = 52 
ts = np.linspace(0,153,n)

# calculate phibar for the different capacity settings
phibar_sigma_plus_low = []
phibar_sigma_plus_high = []
phibar_Omega_low_p_medium = []
phibar_Omega_medium_p_medium = []
phibar_Omega_high_p_medium = []
phibar_Omega_medium_p_low = []
phibar_Omega_medium_p_high = []
phibar_Omega_low_p_high = []
phibar_Omega_high_p_low = []

for o in range(0,n):
    phibar_sigma_plus_low.append(phibar(ts[o],40000,2,150000))
    phibar_sigma_plus_high.append(phibar(ts[o],40000,2,250000))
    phibar_Omega_low_p_medium.append(phibar(ts[o],20000,2,200000))
    phibar_Omega_medium_p_medium.append(phibar(ts[o],40000,2,200000))
    phibar_Omega_high_p_medium.append(phibar(ts[o],60000,2,200000))
    phibar_Omega_medium_p_low.append(phibar(ts[o],40000,1,200000))
    phibar_Omega_medium_p_high.append(phibar(ts[o],40000,10,200000))
    phibar_Omega_low_p_high.append(phibar(ts[o],20000,10,200000))
    phibar_Omega_high_p_low.append(phibar(ts[o],60000,1,200000))

##############
### Figure ###
##############

fig = plt.figure(figsize=(8, 6.5),constrained_layout=True)
grid = fig.add_gridspec(ncols=2, nrows=2, wspace=0.1)
ax1 = fig.add_subplot(grid[0])
ax2 = fig.add_subplot(grid[1],sharey=ax1)
ax3 = fig.add_subplot(grid[3],sharey=ax1)
ax7 = fig.add_subplot(grid[2],sharey=ax1)

ax1.plot(ts,phibar_Omega_low_p_medium,'-',color=My_Blue, linewidth = 2.5)
ax1.plot(ts,phibar_Omega_medium_p_medium,'--',color=My_Amber, linewidth = 2.5)
ax1.plot(ts,phibar_Omega_high_p_medium,'.',color=My_Violet)
ax1.set_xlabel(r'day of intervention $t^*$', fontsize=14)
ax1.set_ylabel('$\overline{\phi}$: maximum allowed level \n of effective contacts',\
               fontsize=14)
ax1.set_xlim(0,153)
ax1.tick_params(axis="x", labelsize=14)
ax1.tick_params(axis="y", labelsize=14)
ax1.legend([r'$\Omega=20000$',r'$\Omega=40000$',r'$\Omega=60000$'],\
           title=r'tracing capacity',\
           prop=dict(size=10.5,weight='bold'))
ax1.text(-0.1, 1.15, string.ascii_uppercase[0], transform=ax1.transAxes, 
        size=20, weight='bold')

ax2.plot(ts,phibar_Omega_medium_p_low,'-',color=My_Blue, linewidth = 2.5)
ax2.plot(ts,phibar_Omega_medium_p_medium,'--',color=My_Amber, linewidth = 2.5)
ax2.plot(ts,phibar_Omega_medium_p_high,'.',color=My_Violet)
ax2.set_xlabel(r'day of intervention $t^*$', fontsize=14)
ax2.set_xlim(0,153)
ax2.tick_params(axis="x", labelsize=14)
ax2.tick_params(axis="y", labelsize=14)
ax2.legend([r'$p=1$',r'$p=2$',r'$p=10$'], title=r'tracing eff. constant',\
           prop=dict(size=10.5,weight='bold'))
ax2.text(-0.1, 1.15, string.ascii_uppercase[1], transform=ax2.transAxes, 
        size=20, weight='bold')

ax3.plot(ts,phibar_Omega_low_p_high,'-',color=My_Blue, linewidth = 2.5)
ax3.plot(ts,phibar_Omega_high_p_low,'--',color=My_Amber, linewidth = 2.5)
ax3.set_xlabel(r'day of intervention $t^*$', fontsize=14)
ax3.set_xlim(0,153)
ax3.tick_params(axis="x", labelsize=14)
ax3.tick_params(axis="y", labelsize=14)
ax3.legend([r'$\Omega=20000,p=10$',r'$\Omega=60000,p=1$'],prop=dict(size=10.5,\
                                                                    weight='bold'))
ax3.text(-0.1, 1.15, string.ascii_uppercase[3], transform=ax3.transAxes, 
        size=20, weight='bold')

ax7.plot(ts,phibar_sigma_plus_low,'-',color=My_Blue, linewidth = 2.5)
ax7.plot(ts,phibar_Omega_medium_p_medium,'--',color=My_Amber, linewidth = 2.5)
ax7.plot(ts,phibar_sigma_plus_high,'.',color=My_Violet)
ax7.set_xlabel(r'day of intervention $t^*$', fontsize=14)
ax7.set_ylabel('$\overline{\phi}$: maximum allowed level \n of effective contacts',\
               fontsize=14)
ax7.set_xlim(0,153)
ax7.tick_params(axis="x", labelsize=14)
ax7.tick_params(axis="y", labelsize=14)
ax7.legend([r'$\sigma_+=150K$',r'$\sigma_+=200K$',r'$\sigma_+=250K$'],\
           title=r'tracing capacity',\
           prop=dict(size=10.5,weight='bold'))
ax7.text(-0.1, 1.15, string.ascii_uppercase[2], transform=ax7.transAxes, 
        size=20, weight='bold')

#plt.savefig('Figure13.pdf', format='pdf', dpi=600)
#plt.savefig('Figure13.eps', format='eps', dpi=600)