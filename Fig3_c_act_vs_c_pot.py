# Figure3: Illustration of c_act as function of c_pot assuming different values 
# for the tracing efficiency constant p

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# costum colors for plotting
My_Amber = "#E69F00"
My_Blue = "#56B4E9"
My_Green = "#009E73"
My_Ver = '#D55E00'

color_list = [My_Blue, My_Amber, My_Green, My_Ver]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)

Omega = 40000 # tracing capacity
c_pot_max = 150000 
c_pot = np.linspace(0,c_pot_max,5*c_pot_max+1);
p_inf = []
p_10 = []
p_2 = []
p_1 = []

for i in range(0,len(c_pot)):
    p_inf.append(c_pot[i]*Omega/np.linalg.norm(np.array([c_pot[i],Omega]),\
                                               ord=np.inf))
    p_10.append(c_pot[i]*Omega/np.linalg.norm(np.array([c_pot[i],Omega]),\
                                              ord=10))
    p_2.append(c_pot[i]*Omega/np.linalg.norm(np.array([c_pot[i],Omega]),\
                                             ord=2))
    p_1.append(c_pot[i]*Omega/np.linalg.norm(np.array([c_pot[i],Omega]),\
                                             ord=1))

##############
### Figure ###
##############

fig = plt.figure(figsize=(4.7, 4),constrained_layout=True)
grid = fig.add_gridspec(ncols=1, nrows=1, wspace=0.1)
ax = fig.add_subplot(grid[0])
ax.plot(c_pot,p_inf,linewidth=2.5)
ax.plot(c_pot,p_10,'--',linewidth=2.5)
ax.plot(c_pot[::1024*20],p_2[::1024*20],'.')
ax.plot(c_pot[::1024*20],p_1[::1024*20],'v',markersize=3.5)
ticks_locx = ax.get_xticks().tolist()
ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_locx)) 
xlabels = ['{:,.0f}'.format(x) + 'K' for x in ax.get_xticks()/1000]
ax.set_xticklabels(xlabels)
ticks_locy = ax.get_yticks().tolist()  
ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_locy)) 
ylabels = ['{:,.0f}'.format(y) + 'K' for y in ax.get_yticks()/1000]
ax.set_yticklabels(ylabels)  
ax.tick_params(axis="x", labelsize=14)
ax.tick_params(axis="y", labelsize=14)
ax.set_ylabel('$c_{\mathrm{act}}$: rate at which contacts \n are isolated',\
              fontsize=14)
ax.set_xlabel('$c_{\mathrm{pot}}$: rate at which contacts \n become traceable',\
              fontsize=14)
ax.set_xlim(0,c_pot_max)
ax.set_ylim(0,Omega+2500)
ax.legend(['p$\\to\infty$','p=10','p=2','p=1'], fontsize=12)

#plt.savefig('Figure3.eps', format='eps', dpi=600)
#plt.savefig('Figure3.pdf', format='pdf', dpi=600)