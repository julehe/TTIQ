""" 
Figure11: Plots showing outcome of simulation of the DDE model in baseline 
setting
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
My_Violet = '#661D98'

color_list = [My_Blue, My_Amber, My_Green, My_Violet]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)

# run simulation 
steps = 32
duration = 130 
T = np.linspace(0,duration,steps*duration+1)
IC = [82975287, 2564, 131, 1301, 86, 52, 2173, 207, 1666, 16533, 0]
def history(t):
    return IC
para_cap = [40000, 2, 200000]
para = [0.6, 0.65, 93, 300, 0.1, 0.2, 2]
para_hist = para
para_interv = para
sol,treff,test_pos,det_ratio,_ \
= model(steps, duration, IC, history, para_cap, para, para_hist, para_interv,\
        0, 1)
# unpack solution
s = np.array(sol[:,0])
e = np.array(sol[:,1])
qe = np.array(sol[:,2])
u1 = np.array(sol[:,3])
qu1 = np.array(sol[:,4])
i1 = np.array(sol[:,5])
u2 = np.array(sol[:,6])
qu2 = np.array(sol[:,7])
i2 = np.array(sol[:,8])
icum = np.array(sol[:,10])
# calcualte total infected and under-ascertainment
u = u1 + u2
qu = qu1 + qu2
i = i1 + i2
infected = e + qe + u + qu + i
under_ascert = (u+qu+i)/(qu+i)
# calculate daily new cases
inewd = np.diff(icum[::steps])/np.diff(T[::steps])
inewd = np.insert(inewd, 0, inewd[0])

##############
### Figure ###
##############

fig = plt.figure(figsize=(9.5, 6.5),constrained_layout=True)
grid = fig.add_gridspec(ncols=2, nrows=2, wspace=0.1)
ax1 = fig.add_subplot(grid[0])
ax2 = fig.add_subplot(grid[1])
ax3 = fig.add_subplot(grid[2])
ax4 = fig.add_subplot(grid[3])

# prevalence
ax1.plot(T,infected,linewidth=2.5)
ax1.set_xlabel('time (days)', fontsize=14)
ax1.set_ylabel('total infected individ. \n (solid)', fontsize=14)
ticks_loc = ax1.get_yticks().tolist()
ax1.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
ylabels = ['{:,.0f}'.format(y) + 'K' for y in ax1.get_yticks()/1000]
ax1.set_yticklabels(ylabels)   
ax1.tick_params(axis="x", labelsize=14)
ax1.tick_params(axis="y", labelsize=14)
ax1.set_xlim(0,duration)
ax1.text(-0.1, 1.15, string.ascii_uppercase[0], transform=ax1.transAxes, 
        size=20, weight='bold')
# susceptibles
ax12=ax1.twinx()
ax12.plot(T,s,'--',linewidth=2.5,color=My_Amber)
ax12.set_xlabel('time (days)', fontsize=14)
ax12.set_ylabel('susceptibles \n (dashed)', fontsize=14)
ax12.tick_params(axis="x", labelsize=14)
ax12.tick_params(axis="y", labelsize=14)

# daily incidence
ax2.bar(T[::steps],inewd,width=1)
ax2.set_xlabel('time (days) $t$', fontsize=14)
ax2.set_ylabel('daily case incidence \n (solid)', fontsize=14)  
ticks_loc = ax2.get_yticks().tolist()
ax2.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc)) 
ylabels = ['{:,.0f}'.format(y) + 'K' for y in ax2.get_yticks()/1000]
ax2.set_yticklabels(ylabels)  
ax2.tick_params(axis="x", labelsize=14)
ax2.tick_params(axis="y", labelsize=14)
ax2.set_xlim(0,duration)
ax2.text(-0.1, 1.15, string.ascii_uppercase[1], transform=ax2.transAxes, 
        size=20, weight='bold')

#tracing efficiency
ax22=ax2.twinx()
ax22.plot(T,treff,'--',linewidth=2.5,color=My_Amber)
ax22.set_xlabel('time (days) $t$', fontsize=14)
ax22.set_ylabel('tracing efficiency \n (dashed)', fontsize=14)
ax22.tick_params(axis="x", labelsize=14)
ax22.tick_params(axis="y", labelsize=14)

# test positive rate
ax3.plot(T,test_pos,linewidth=2.5)
ax3.set_xlabel('time (days)', fontsize=14)
ax3.set_ylabel('test positive rate \n (solid)', fontsize=14)
ax3.tick_params(axis="x", labelsize=14)
ax3.tick_params(axis="y", labelsize=14)    
ax3.set_xlim(0,duration)
ax3.text(-0.1, 1.15, string.ascii_uppercase[2], transform=ax3.transAxes, 
        size=20, weight='bold')

# detection ratio
ax32=ax3.twinx()
ax32.plot(T,det_ratio,'--',linewidth=2.5,color=My_Amber)
ax32.set_xlabel('time (days) $t$', fontsize=14)
ax32.set_ylabel('detection ratio \n (dashed)', fontsize=14)
ax32.tick_params(axis="x", labelsize=14)
ax32.tick_params(axis="y", labelsize=14)


# under-ascertainment
ax4.plot(T,(u+qu+i)/(qu+i),linewidth=2.5)
ax4.set_xlabel('time (days)', fontsize=14)#
ax4.set_ylabel('underascertainment', fontsize=14)
ax4.tick_params(axis="x", labelsize=14)
ax4.tick_params(axis="y", labelsize=14)
ax4.set_xlim(0,duration)
ax4.text(-0.1, 1.15, string.ascii_uppercase[3], transform=ax4.transAxes, 
        size=20, weight='bold')
ax42=ax4.twinx()
ax42.set_yticks([])
fig.align_ylabels()
plt.show()

#fig.savefig('Figure11.eps', format='eps', dpi=600)
#fig.savefig('Figure11.pdf', format='pdf', dpi=600)
