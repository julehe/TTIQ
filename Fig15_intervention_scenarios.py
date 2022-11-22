"""
Figure 14: simulation of the DDE model with different interventions
"""

import numpy as np
import matplotlib.pyplot as plt
import string
import matplotlib.ticker as mticker
from model import model

# costum colors for plotting
My_Amber = "#E69F00"
My_Blue = "#56B4E9"
My_Green = "#009E73"
My_Ver = '#D55E00'

color_list = [My_Blue, My_Amber, My_Green, My_Ver]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)

# tick formatter: Turns large tick values 
# (in the billions, millions and thousands) into xK, xM, xB format
def reformat_large_tick_values(tick_val, pos):
    if tick_val >= 1000000000:
        val = round(tick_val/1000000000, 1)
        new_tick_format = '{:}B'.format(val)
    elif tick_val >= 1000000:
        val = round(tick_val/1000000, 1)
        new_tick_format = '{:}M'.format(val)
    elif tick_val >= 1000:
        val = round(tick_val/1000, 1)
        new_tick_format = '{:}K'.format(val)
    elif tick_val < 1000:
        new_tick_format = round(tick_val, 1)
    else:
        new_tick_format = tick_val
    # make new_tick_format into a string value
    new_tick_format = str(new_tick_format)
    # code below will keep 4.5M as is but change values such as 4.0M to 4M 
    # since that zero after the decimal isn't needed
    index_of_decimal = new_tick_format.find(".")
    
    if index_of_decimal != -1:
        value_after_decimal = new_tick_format[index_of_decimal+1]
        if value_after_decimal == "0":
            # remove the 0 after the decimal point since it's not needed
            new_tick_format = new_tick_format[0:index_of_decimal]\
                              + new_tick_format[index_of_decimal+2:]
    return new_tick_format

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

steps = 32;
duration = 188
T = np.linspace(0,duration,steps*duration+1);
IC = [82975287, 2564, 131, 1301, 86, 52, 2173, 207, 1666, 16533, 0]
def history(t):
    return IC

def interv_simul(interv_time, trcov, sigma_U, phi):
    para = [phi_b, trcov_b, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b,tau_b]
    para_hist = para
    para_interv = [phi, trcov, sigma_U, sigma_Q_b, p_I_b, p_Q_b,tau_b]
    sol,*_ =\
    model(steps, duration, IC, history, [40000,2,200000], para, para_hist,\
          para_interv, interv_time, 0)
    e = np.array(sol[:,1]);
    qe = np.array(sol[:,2]);
    u1 = np.array(sol[:,3]);
    u2 = np.array(sol[:,6]);
    qu1 = np.array(sol[:,4]);
    qu2 = np.array(sol[:,7]);
    i1 = np.array(sol[:,5]);
    i2 = np.array(sol[:,8]);
    inf = e + qe + u1 + u2 + qu1 + qu2 + i1 + i2
    return inf

# run simulation with different interventions and different timing
early = 47
late = 123

phi_bas = 0.49
phi_strict = 0.46

# only reducing contacts 
inf_only_phi_early = interv_simul(early, trcov_b, sigma_U_b, phi_bas) # early
inf_only_phi_late = interv_simul(late, trcov_b, sigma_U_b, phi_bas) # late

# reducing contacts and improved tracing coverage
inf_with_trcov_early = interv_simul(early, 1, sigma_U_b, phi_bas) # early
inf_with_trcov_late = interv_simul(late, 1, sigma_U_b, phi_bas) # late

# reducing contacts and improved testing
inf_with_sigmaU_early = interv_simul(early, trcov_b, 118, phi_bas) # early
inf_with_sigmaU_late = interv_simul(late, trcov_b, 118, phi_bas) # late

# stricter reduction of contacts
inf_strict_phi_early = interv_simul(early, trcov_b, sigma_U_b, phi_strict) # early
inf_strict_phi_late = interv_simul(late, trcov_b, sigma_U_b, phi_strict) # late

##############
### FIGURE ###
##############

plt.rcParams['legend.title_fontsize'] = 'large'

fig = plt.figure(figsize=(9.7, 4),constrained_layout=True)
grid = fig.add_gridspec(ncols=2, nrows=1, wspace=0.1)
ax1 = fig.add_subplot(grid[0])
ax2 = fig.add_subplot(grid[1])

ax1.plot(T,inf_only_phi_early,linewidth=2.5)
ax1.plot(T,inf_with_trcov_early,'--',linewidth=2.5)
ax1.plot(T[::5*steps],inf_with_sigmaU_early[::5*steps],'.',linewidth=2.5)
ax1.plot(T[::7*steps],inf_strict_phi_early[::7*steps],'v',markersize=3.5)
ax1.tick_params(axis="x", labelsize=14)
ax1.tick_params(axis="y", labelsize=14)
ax1.set_xlabel(r'time (days)', fontsize=14)
ax1.set_ylabel('total infected individuals', fontsize=14)
ax1.yaxis.set_major_formatter(mticker.FuncFormatter(reformat_large_tick_values)); 
ax1.set_xlim(0,duration)
ax1.set_title('early intervention', fontsize=16)
ax1.legend([r'$\phi_1=0.49$',r'$\phi_1=0.49,\omega_1=1.0$',\
            r'$\phi_1=0.49,\sigma_{U_2,1}=118$',r'$\phi_1=0.46$'],\
            prop=dict(size=10.5,weight='bold'))
ax1.text(-0.1, 1.15, string.ascii_uppercase[0], transform=ax1.transAxes, 
        size=20, weight='bold')


ax2.plot(T,inf_only_phi_late,linewidth=2.5)
ax2.plot(T,inf_with_trcov_late,'--',linewidth=2.5)
ax2.plot(T[::5*steps],inf_with_sigmaU_late[::5*steps],'.',linewidth=2.5)
ax2.plot(T[::7*steps],inf_strict_phi_late[::7*steps],'v',markersize=3.5)
ax2.tick_params(axis="x", labelsize=14)
ax2.tick_params(axis="y", labelsize=14)
ax2.set_xlabel(r'time (days)', fontsize=14)
ax2.yaxis.set_major_formatter(mticker.FuncFormatter(reformat_large_tick_values)); 
ax2.set_xlim(0,duration)
ax2.set_title('late intervention', fontsize=16)
ax2.legend([r'$\phi_1=0.49$',r'$\phi_1=0.49,\omega_1=1.0$',\
            r'$\phi_1=0.49,\sigma_{U_2,1}=118$',r'$\phi_1=0.46$'],\
            prop=dict(size=10.5,weight='bold'))
ax2.text(-0.1, 1.15, string.ascii_uppercase[1], transform=ax2.transAxes, 
        size=20, weight='bold')
fig.align_ylabels()

#fig.savefig('Figure14.eps', format='eps', dpi=600)
#fig.savefig('Figure14.pdf', format='pdf', dpi=600)

