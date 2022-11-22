"""
Figure 13: $\bar\phi$ at an early and a late intervention time point + varying
individual TTIQ parameters  
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import string
from model import model

# costum colors for plotting 
My_Amber = "#E69F00"
My_Blue = "#56B4E9"
My_Green = "#009E73"
My_Violet = '#661D98'

color_list = [My_Blue, My_Amber, My_Green, My_Violet]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)

#baseline parameter setting and considered intervals
# level of effective contacts
phi_b = 0.6
# tracing coverage
trcov_vec = np.linspace(0,1,20)
trcov_b = 0.65
# rel. freq. of testing in U
sigma_U_vec = np.linspace(1,186,20)
sigma_U_b = 93
# rel. freq. of testing in Q
sigma_Q_vec = np.linspace(1,600,35)
sigma_Q_b = 300
# strictness of isolation 
p_I_vec = np.linspace(0,1,20)
p_I_b = 0.1
# strictness of quarantine
p_Q_vec = np.linspace(0,1,20)
p_Q_b = 0.2
# tracing delay
tau_vec = np.linspace(0.5,14,20)
tau_b = 2
# first running the simulation in baseline setting on total time horizon.
# saves computing time since later simulations can use interpolated solution
# as history and do not have to run whole simulation again.
steps = 32
duration = 200
T = np.linspace(0,duration,steps*duration+1);
IC = [82975287, 2564, 131, 1301, 86, 52, 2173, 207, 1666, 16533, 0]
def history(t):
    return IC
para_cap = [40000, 2, 200000]
para = [phi_b, trcov_b, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b,tau_b]
para_hist = para
para_interv = para
sol,*_ =\
model(steps, duration, IC, history, para_cap, para, para_hist, para_interv, 0, 0)

s = interp1d(T,np.array(sol[:,0]), kind='linear')
e = interp1d(T,np.array(sol[:,1]), kind='linear')
qe = interp1d(T,np.array(sol[:,2]), kind='linear')
u1 = interp1d(T,np.array(sol[:,3]), kind='linear')
qu1 = interp1d(T,np.array(sol[:,4]), kind='linear')
i1 = interp1d(T,np.array(sol[:,5]), kind='linear')
u2 = interp1d(T,np.array(sol[:,6]), kind='linear')
qu2 = interp1d(T,np.array(sol[:,7]), kind='linear')
i2 = interp1d(T,np.array(sol[:,8]), kind='linear')
r = interp1d(T,np.array(sol[:,9]), kind='linear')

# phibar for a given intervention timing and TTIQ paramater setting
def phibar(trcov, sigma_U, sigma_Q, p_I, p_Q, tau, timing):
    upper = 1
    lower = 0
    middle = lower+(upper-lower)/2
    while True:
        steps = 32
        duration_temp = 56
        phi = middle
        para = [phi, trcov, sigma_U, sigma_Q, p_I, p_Q, tau]
        para_hist = [phi_b, trcov_b, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b,tau_b]
        para_interv = para
        #history function
        def history(t):
            if timing == "early":
                t = 47 + t
            else: 
                t = 123 + t
            return [s(t),e(t),qe(t),u1(t),qu1(t),i1(t),u2(t),qu2(t),i2(t),r(t),0]
        *_,lamS =\
        model(steps,duration_temp,IC,history,para_cap,para,para_hist,para_interv,0,0)
        if (lamS)[14*steps]>=max((lamS)[14*steps:]) or (lamS)[0]>=max((lamS)[0:]):    
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

# calculate phibar in the different settings
phibar_trcov_early = []

phibar_sigma_U_early = []

phibar_sigma_Q_early = []

phibar_p_I_early = []

phibar_p_Q_early = []

phibar_tau_early = []

for i in trcov_vec:
    phibar_trcov_early.append(phibar(i, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b,\
                                     tau_b, 'early'))
   

for i in sigma_U_vec:
    phibar_sigma_U_early.append(phibar(trcov_b, i, sigma_Q_b, p_I_b, p_Q_b,\
                                       tau_b,'early'))

for i in sigma_Q_vec:
    phibar_sigma_Q_early.append(phibar(trcov_b, sigma_U_b, i, p_I_b, p_Q_b,\
                                       tau_b,'early'))
    
for i in p_I_vec:
    phibar_p_I_early.append(phibar(trcov_b, sigma_U_b, sigma_Q_b, i, p_Q_b,\
                                   tau_b,'early'))
  
for i in p_Q_vec:
    phibar_p_Q_early.append(phibar(trcov_b, sigma_U_b, sigma_Q_b, p_I_b, i,\
                                   tau_b,'early'))

for i in tau_vec:
    phibar_tau_early.append(phibar(trcov_b, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b,\
                                   i,'early'))

phibar_trcov_late = []

phibar_sigma_U_late = []

phibar_sigma_Q_late = []

phibar_p_I_late = []

phibar_p_Q_late = []

phibar_tau_late = []

for i in trcov_vec:
    phibar_trcov_late.append(phibar(i, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b,\
                                    tau_b, 'late'))
    

for i in sigma_U_vec:
    phibar_sigma_U_late.append(phibar(trcov_b, i, sigma_Q_b, p_I_b, p_Q_b,\
                                      tau_b,'late'))

for i in sigma_Q_vec:
    phibar_sigma_Q_late.append(phibar(trcov_b, sigma_U_b, i, p_I_b, p_Q_b,\
                                      tau_b,'late'))
    
for i in p_I_vec:
    phibar_p_I_late.append(phibar(trcov_b, sigma_U_b, sigma_Q_b, i, p_Q_b,\
                                  tau_b,'late'))

for i in p_Q_vec:
    phibar_p_Q_late.append(phibar(trcov_b, sigma_U_b, sigma_Q_b, p_I_b, i,\
                                  tau_b,'late'))
  
for i in tau_vec:
    phibar_tau_late.append(phibar(trcov_b, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b,\
                                  i,'late'))
   
##############
### Figure ###
##############

fig = plt.figure(figsize=(14, 8.5),constrained_layout=True)
grid = fig.add_gridspec(ncols=3, nrows=2, wspace=0.1)
ax1 = fig.add_subplot(grid[0])
ax2 = fig.add_subplot(grid[1])
ax3 = fig.add_subplot(grid[2])
ax4 = fig.add_subplot(grid[3])
ax5 = fig.add_subplot(grid[4])
ax6 = fig.add_subplot(grid[5])


ax1.plot(p_I_vec,phibar_p_I_early,'-',p_I_vec,phibar_p_I_late,'--',linewidth=2.5)
ax1.set_xlabel(r'$p_I$: strictness of isolation in $I$', fontsize=14)
ax1.set_xlim(np.min(p_I_vec),np.max(p_I_vec))
ax1.set_ylim(0.3,0.6)
ax1.tick_params(axis="x", labelsize=14)
ax1.tick_params(axis="y", labelsize=14)
ax1.legend(['early intervention','late intervention'],prop=dict(size=10.5,\
                                                                weight='bold'))
ax1.text(-0.1, 1.15, string.ascii_uppercase[0], transform=ax1.transAxes, 
        size=20, weight='bold')

ax2.plot(sigma_U_vec,phibar_sigma_U_early,'-',sigma_U_vec,phibar_sigma_U_late,\
         '--',linewidth=2.5)
ax2.set_xlabel(r'$\sigma_{U_2}$: relative frequency of testing in $U_2$', fontsize=14)
ax2.set_xlim(np.min(sigma_U_vec),np.max(sigma_U_vec))
ax2.set_ylim(0.3,0.6)
ax2.tick_params(axis="x", labelsize=14)
ax2.tick_params(axis="y", labelsize=14)
ax2.legend(['early intervention','late intervention'],prop=dict(size=10.5,\
                                                                weight='bold'))
ax2.text(-0.1, 1.15, string.ascii_uppercase[1], transform=ax2.transAxes, 
        size=20, weight='bold')

ax3.plot(trcov_vec,phibar_trcov_early,'-',trcov_vec,phibar_trcov_late,'--',\
         linewidth=2.5)
ax3.set_xlabel(r'$\omega$: tracing coverage', fontsize=14)
ax3.set_xlim(0,1)
ax3.set_ylim(0.3,0.6)
ax3.tick_params(axis="x", labelsize=14)
ax3.tick_params(axis="y", labelsize=14)
ax3.legend(['early intervention','late intervention'],prop=dict(size=10.5,\
                                                                weight='bold'))
ax3.text(-0.1, 1.15, string.ascii_uppercase[2], transform=ax3.transAxes, 
        size=20, weight='bold')

ax4.plot(tau_vec,phibar_tau_early,'-',tau_vec,phibar_tau_late,'--',linewidth=2.5)
ax4.set_xlabel(r'$\kappa$: tracing delay', fontsize=14)
ax4.set_xlim(np.min(tau_vec),np.max(tau_vec))
ax4.set_ylim(0.3,0.6)
ax4.tick_params(axis="x", labelsize=14)
ax4.tick_params(axis="y", labelsize=14)
ax4.legend(['early intervention','late intervention'],prop=dict(size=10.5,\
                                                                weight='bold'))
ax4.text(-0.1, 1.15, string.ascii_uppercase[3], transform=ax4.transAxes, 
        size=20, weight='bold')
fig.align_ylabels()

ax5.plot(sigma_Q_vec,phibar_sigma_Q_early,'-',sigma_Q_vec,phibar_sigma_Q_late,\
         '--',linewidth=2.5)
ax5.set_xlabel(r'$\sigma_{Q}$: relative frequency of testing in $Q$', fontsize=14)
ax5.set_xlim(np.min(sigma_Q_vec),np.max(sigma_Q_vec))
ax5.set_ylim(0.3,0.6)
ax5.tick_params(axis="x", labelsize=14)
ax5.tick_params(axis="y", labelsize=14)
ax5.legend(['early intervention','late intervention'],prop=dict(size=10.5,\
                                                                weight='bold'))
ax5.text(-0.1, 1.15, string.ascii_uppercase[4], transform=ax5.transAxes, 
        size=20, weight='bold')

ax6.plot(p_Q_vec,phibar_p_Q_early,'-',p_Q_vec,phibar_p_Q_late,'--',linewidth=2.5)
ax6.set_xlabel(r'$p_Q$: strictness of isolation in $Q$', fontsize=14)
ax6.set_xlim(np.min(p_Q_vec),np.max(p_Q_vec))
ax6.set_ylim(0.3,0.6)
ax6.tick_params(axis="x", labelsize=14)
ax6.tick_params(axis="y", labelsize=14)
ax6.legend(['early intervention','late intervention'],prop=dict(size=10.5,\
                                                                weight='bold'))
ax6.text(-0.1, 1.15, string.ascii_uppercase[5], transform=ax6.transAxes, 
        size=20, weight='bold')
#plt.savefig('Figure13.eps',format='eps',dpi=600)
#plt.savefig('Figure13.pdf',format='pdf',dpi=600)



