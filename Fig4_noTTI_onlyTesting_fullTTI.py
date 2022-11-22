"""
Figure 4: Comparison of $\phi^*$ in scenario without any TTIQ (no TTIQ) vs.
only testing (only testing) vs. testing and tracing (full TTIQ).  
+ The same with improved testing.
"""

import matplotlib.pyplot as plt
import string
from phistar import phistar

# costum colors for plotting 
My_Amber = "#E69F00"
My_Blue = "#56B4E9"
My_Green = "#009E73"
My_Ver = '#D55E00'

color_list = [My_Amber, My_Blue, My_Green, My_Ver]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)

# baseline parameter setting
cov_b = 0.65 # tracing coverage
sigma_U_b = 93 # rel. freq. of testing in U2
sigma_Q_b = 300 # rel. freq. of testing in Q
p_I_b = 0.1 # strictness of isolation in I
p_Q_b = 0.2 # strictness of quarantine in Q
tau_b = 2 # tracing delay
theta_b = 1.5 # scaling factor of presymptomatic transmission rate
lat_b = 3.5 # average duraiton of latency period
pre_b = 2 # average duration of presymptomatic infectious period
inf_b = 9 # average total duration of infectious period
R0 = 3.3 # basic reproduction number

# to consider the scneario with no TTI we set p_Q and p_I to 1
# to consider the scenario with only testing but no tracing we set the tracing
# coverage to 0

no_TTI = phistar([cov_b,sigma_U_b,sigma_Q_b,1,1,tau_b,theta_b,lat_b,pre_b,\
                  inf_b,R0],0)
only_testing = phistar([0,sigma_U_b,sigma_Q_b,p_I_b,p_Q_b,tau_b,theta_b,lat_b,\
                        pre_b,inf_b,R0],0)
full_TTI = phistar([cov_b,sigma_U_b,sigma_Q_b,p_I_b,p_Q_b,tau_b,theta_b,lat_b,\
                    pre_b,inf_b,R0],0)
phistars = [no_TTI, only_testing, full_TTI]

# now the same with improved testing 

sigma_U_improv = 185

no_TTI_improv = phistar([cov_b,sigma_U_improv,sigma_Q_b,1,1,tau_b,theta_b,\
                         lat_b,pre_b,inf_b,R0],0)
only_testing_improv = phistar([0,sigma_U_improv,sigma_Q_b,p_I_b,p_Q_b,tau_b,\
                               theta_b,lat_b,pre_b,inf_b,R0],0)
full_TTI_improv = phistar([cov_b,sigma_U_improv,sigma_Q_b,p_I_b,p_Q_b,tau_b,\
                           theta_b,lat_b,pre_b,inf_b,R0],0)
phistars_improv = [no_TTI_improv, only_testing_improv, full_TTI_improv]

############
## FIGURE ##
############
# labels for x-axis
names = ["no TTIQ", "only testing", "full TTIQ"]

fig = plt.figure(figsize=(9.8, 4.5),constrained_layout=True)
grid = fig.add_gridspec(ncols=2, nrows=1, wspace=0.1)
ax1 = fig.add_subplot(grid[0])
ax2 = fig.add_subplot(grid[1])

ax1.bar(names,phistars,width=0.4,color=[My_Blue,My_Blue,My_Blue])
ax1.set_ylabel('$\\phi^*$: maximum allowed level of \n effective contacts',\
               fontsize=14)  
ax1.set_ylim(0,0.61)
ax1.tick_params(axis="x", labelsize=14)
ax1.tick_params(axis="y", labelsize=14)
ax1.set_title('baseline setting', fontsize=16)
ax1.annotate(str(round(phistars[0],3)), xy=(names[0],phistars[0]),\
             ha='center', va='bottom', fontsize=14)
ax1.annotate(str(round(phistars[1],3)), xy=(names[1],phistars[1]),\
             ha='center', va='bottom', fontsize=14)
ax1.annotate(str(round(phistars[2],3)), xy=(names[2],phistars[2]),\
             ha='center', va='bottom', fontsize=14)
ax1.text(-0.1, 1.15, string.ascii_uppercase[0], transform=ax1.transAxes, 
        size=20, weight='bold')


ax2.bar(names,phistars_improv,width=0.4,color=[My_Blue,My_Blue,My_Blue])
ax2.set_ylim(0,0.61)
ax2.tick_params(axis="x", labelsize=14)
ax2.tick_params(axis="y", labelsize=14)
ax2.set_title('setting with improved testing', fontsize=16)
ax2.annotate(str(round(phistars_improv[0],3)),\
             xy=(names[0],phistars_improv[0]), ha='center',va='bottom',\
             fontsize=14)
ax2.annotate(str(round(phistars_improv[1],3)),\
             xy=(names[1],phistars_improv[1]), ha='center', va='bottom',\
            fontsize=14)
ax2.annotate(str(round(phistars_improv[2],3)),\
             xy=(names[2],phistars_improv[2]), ha='center', va='bottom',\
             fontsize=14)
ax2.text(-0.1, 1.15, string.ascii_uppercase[1], transform=ax2.transAxes, 
        size=20, weight='bold')
fig.align_ylabels()

#plt.savefig('Figure5.eps', format='eps', dpi=600)
#plt.savefig('Figure5.pdf', format='pdf', dpi=600)




