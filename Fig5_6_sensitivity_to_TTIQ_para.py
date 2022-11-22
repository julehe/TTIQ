"""
Figure 6 & 7: phistar varying TTIQ parameters.
"""

import numpy as np
import matplotlib.pyplot as plt
from phistar import phistar
import string


# costum colors for plotting
My_Amber = "#E69F00"
My_Blue = "#56B4E9"
My_Green = "#009E73"
My_Violet = '#661D98'

color_list = [My_Blue, My_Amber, My_Green, My_Violet]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)

# baseline parameter setting + ranges considered for sensitivity wrt TTIQ paras
cov_b = 0.65 # tracing coverage
cov = np.linspace(0,1,51)
sigma_U_b = 93 # rel. freq. of testing in U2
sigma_U = np.linspace(1,186,93)
sigma_Q_b = 300 # rel. freq. of testing in Q
sigma_Q = np.linspace(1,600,150)
p_I_b = 0.1 # strictness of isolation
p_I = np.linspace(0,1,51)
p_Q_b = 0.2 # strictness of quarantine
p_Q = np.linspace(0,1,51)
tau_b = 2 # tracing delay
tau = np.linspace(0.5,14,50)
theta_b = 1.5 # scaling factor of presymptomatic transmission rate
lat_b = 3.5 # average duraiton of latency period
pre_b = 2 # average duration of presymptomatic infectious period
inf_b = 9 # average total duration of infectious period
R0 = 3.3 # basic reproduction number

phistar_cov = []

phistar_sigma_U = []

phistar_sigma_Q = []

phistar_p_Q = []

phistar_p_I = []

phistar_tau = []

phistar_baseline = phistar([cov_b, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b, tau_b,\
                            theta_b, lat_b, pre_b, inf_b, R0],0)

##############################################
### varying TTIQ paras one after the other ###
##############################################

for i in cov:
    phistar_cov.append(phistar([i, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b, tau_b,\
                                theta_b, lat_b, pre_b, inf_b, R0],0))
   

for i in sigma_U:
    phistar_sigma_U.append(phistar([cov_b, i, sigma_Q_b, p_I_b, p_Q_b, tau_b,\
                                    theta_b, lat_b, pre_b, inf_b, R0],0))
   
for i in sigma_Q:
    phistar_sigma_Q.append(phistar([cov_b, sigma_U_b, i, p_I_b, p_Q_b, tau_b,\
                                    theta_b, lat_b, pre_b, inf_b, R0],0))
   
for i in p_Q:
    phistar_p_Q.append(phistar([cov_b, sigma_U_b, sigma_Q_b, p_I_b, i, tau_b,\
                                theta_b, lat_b, pre_b, inf_b, R0],0))
  
for i in p_I:
    phistar_p_I.append(phistar([cov_b, sigma_U_b, sigma_Q_b, i, p_Q_b, tau_b,\
                                theta_b, lat_b, pre_b, inf_b, R0],0))
  
for i in tau:
    phistar_tau.append(phistar([cov_b, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b, i,\
                                theta_b, lat_b, pre_b, inf_b, R0],0))

#####################
### p_Q & sigma_Q ###
#####################

n=100
p_Q_heat, sigma_Q_heat = np.meshgrid(np.linspace(0.1, 1, n), np.linspace(1,600, n))
phistar_p_Q_sigma_Q = []
for o in range(0,n):
    entry = []
    p_Q_heat_temp = p_Q_heat[1]
    sigma_Q_heat_temp = sigma_Q_heat[:,1]
    for j in range(0,n):
        entry.append(phistar([cov_b, sigma_U_b, sigma_Q_heat_temp[o], p_I_b,\
                              p_Q_heat_temp[j], tau_b, theta_b, lat_b,\
                              pre_b, inf_b, R0],0))
    phistar_p_Q_sigma_Q .append(entry)
    
phistar_p_Q_sigma_Q = np.array(phistar_p_Q_sigma_Q)
l_p_Q_heat = p_Q_heat.min()
r_p_Q_heat = p_Q_heat.max()
l_sigma_Q_heat = sigma_Q_heat.min()
r_sigma_Q_heat = sigma_Q_heat.max()
l_c_p_Q_sigma_Q, r_c_p_Q_sigma_Q  = np.abs(phistar_p_Q_sigma_Q).min(),\
                                    np.abs(phistar_p_Q_sigma_Q).max()


#############################################################################
### sigma_U in (1) low compliance seeting vs. (2) high compliance setting ###
#############################################################################

phistar_sigma_U_low_compl = []
for i in sigma_U:
    phistar_sigma_U_low_compl.append(phistar([cov_b, i, sigma_Q_b, 0.35, 0.7,\
                                              tau_b, theta_b, lat_b, pre_b,\
                                              inf_b, R0],0)) 
     

#####################
### cov & sigma_U ###
#####################

cov_heat, sigma_U_heat = np.meshgrid(np.linspace(0, 1, n), np.linspace(0,186, n))
phistar_cov_sigma_U = []
for o in range(0,n):
    entry = [];
    cov_heat_temp = cov_heat[1];
    sigma_U_heat_temp = sigma_U_heat[:,1];
    for j in range(0,n):
        entry.append(phistar([cov_heat_temp[j], sigma_U_heat_temp[o],\
                              sigma_Q_b, p_I_b, p_Q_b, tau_b, theta_b, lat_b,\
                              pre_b, inf_b, R0],0));
    phistar_cov_sigma_U.append(entry)
    
phistar_cov_sigma_U = np.array(phistar_cov_sigma_U)
l_cov_heat = cov_heat.min()
r_cov_heat = cov_heat.max()
l_sigma_U_heat = sigma_U_heat.min()
r_sigma_U_heat = sigma_U_heat.max()
l_c_cov_sigma_U,r_c_cov_sigma_U  = np.abs(phistar_cov_sigma_U).min(),\
                                   np.abs(phistar_cov_sigma_U).max()


##############
### FIGURE ###
##############

cmap = plt.get_cmap('viridis')
legend_properties = {'weight':'bold'}

fig = plt.figure(figsize=(14, 8.5),constrained_layout=True)
grid = fig.add_gridspec(ncols=3, nrows=2, wspace=0.1)
ax1 = fig.add_subplot(grid[0])
ax3 = fig.add_subplot(grid[1])
ax4 = fig.add_subplot(grid[2])
ax5 = fig.add_subplot(grid[3])
ax6 = fig.add_subplot(grid[4])
ax7 = fig.add_subplot(grid[5])

ax1.plot(p_I,phistar_p_I,linewidth=2.5)
ax1.plot((0, p_I[-1]), (phistar_baseline, phistar_baseline), '--',\
          color='dimgrey',linewidth=2)
ax1.set_xlabel(r'$p_I$: strictness of isolation in $I$', fontsize=14)
ax1.set_xlim(np.min(p_I),np.max(p_I))
ax1.set_ylim(0.3,0.6)
ax1.tick_params(axis="x", labelsize=14)
ax1.tick_params(axis="y", labelsize=14)
ax1.text(0.7,0.55,'baseline',transform=ax1.transAxes, 
        size=13, weight='bold',color='dimgrey')
ax1.text(-0.1, 1.15, string.ascii_uppercase[0], transform=ax1.transAxes, 
        size=20, weight='bold')

ax3.plot(sigma_U,phistar_sigma_U,linewidth=2.5)
ax3.plot((0, sigma_U[-1]), (phistar_baseline, phistar_baseline), '--',\
          color='dimgrey',linewidth=2)
ax3.set_xlabel(r'$\sigma_{U_2}$: relative frequency of testing in $U_2$',\
               fontsize=14) 
ax3.set_xlim(np.min(sigma_U),np.max(sigma_U))
ax3.tick_params(axis="x", labelsize=14)
ax3.tick_params(axis="y", labelsize=14)
ax3.set_ylim(0.3,0.6)
ax3.text(0.05,0.55,'baseline',transform=ax3.transAxes, 
        size=13, weight='bold',color='dimgrey')
ax3.text(-0.1, 1.15, string.ascii_uppercase[1], transform=ax3.transAxes, 
        size=20, weight='bold')

ax4.plot(cov,phistar_cov,linewidth=2.5)
ax4.plot((0, cov[-1]), (phistar_baseline, phistar_baseline), '--',\
          color='dimgrey',linewidth=2)
ax4.set_xlabel(r'$\omega$: tracing coverage', fontsize=14) 
ax4.set_xlim(0,1)
ax4.tick_params(axis="x", labelsize=14)
ax4.tick_params(axis="y", labelsize=14)
ax4.set_ylim(0.3,0.6)
ax4.text(0.05,0.55,'baseline',transform=ax4.transAxes, 
        size=13, weight='bold',color='dimgrey')
ax4.text(-0.1, 1.15, string.ascii_uppercase[2], transform=ax4.transAxes, 
        size=20, weight='bold',color='dimgrey')

ax5.plot(tau,phistar_tau,linewidth=2.5)
ax5.plot((0, tau[-1]), (phistar_baseline, phistar_baseline), '--',\
          color='dimgrey',linewidth=2)
ax5.set_xlabel(r'$\kappa$: tracing delay', fontsize=14)
ax5.set_xlim(np.min(tau),np.max(tau))
ax5.set_ylim(0.3,0.6)
ax5.tick_params(axis="x", labelsize=14)
ax5.tick_params(axis="y", labelsize=14)
ax5.text(0.7,0.55,'baseline',transform=ax5.transAxes, 
        size=13, weight='bold',color='dimgrey')
ax5.text(-0.1, 1.15, string.ascii_uppercase[3], transform=ax5.transAxes, 
        size=20, weight='bold')

ax6.plot(sigma_Q,phistar_sigma_Q,linewidth=2.5)
ax6.plot((0, sigma_Q[-1]), (phistar_baseline, phistar_baseline), '--',\
          color='dimgrey',linewidth=2)
ax6.set_xlabel(r'$\sigma_{Q}$: relative frequency of testing in $Q$',\
               fontsize=14) 
ax6.set_xlim(np.min(sigma_Q),np.max(sigma_Q))
ax6.set_ylim(0.3,0.6)
ax6.tick_params(axis="x", labelsize=14)
ax6.tick_params(axis="y", labelsize=14)
ax6.text(0.7,0.55,'baseline',transform=ax6.transAxes, 
        size=13, weight='bold',color='dimgrey')
ax6.text(-0.1, 1.15, string.ascii_uppercase[4], transform=ax6.transAxes, 
        size=20, weight='bold')

ax7.plot(p_Q,phistar_p_Q,linewidth=2.5)
ax7.plot((0, p_Q[-1]), (phistar_baseline, phistar_baseline), '--',\
          color='dimgrey',linewidth=2)
ax7.set_xlabel(r'$p_Q$: strictness of isolation in $Q$', fontsize=14)
ax7.set_xlim(np.min(p_Q),np.max(p_Q))
ax7.set_ylim(0.3,0.6)
ax7.tick_params(axis="x", labelsize=14)
ax7.tick_params(axis="y", labelsize=14)
ax7.text(0.7,0.55,'baseline',transform=ax7.transAxes, 
        size=13, weight='bold',color='dimgrey')
ax7.text(-0.1, 1.15, string.ascii_uppercase[5], transform=ax7.transAxes, 
        size=20, weight='bold')



#plt.savefig('Figure6.eps', format='eps', dpi=600)
#plt.savefig('Figure6.pdf', format='pdf', dpi=600)

fig2 = plt.figure(figsize=(16.5, 4.7),constrained_layout=True)
grid = fig2.add_gridspec(ncols=3, nrows=1, wspace=0.1)
ax10 = fig2.add_subplot(grid[0])
ax8 = fig2.add_subplot(grid[1])
ax11 = fig2.add_subplot(grid[2])

pl10 = ax10.pcolormesh(p_Q_heat, sigma_Q_heat, phistar_p_Q_sigma_Q, cmap=cmap,\
                       vmin=l_c_p_Q_sigma_Q, vmax=r_c_p_Q_sigma_Q, shading='auto')
ax10.axis([l_p_Q_heat, r_p_Q_heat, l_sigma_Q_heat, r_sigma_Q_heat])
ax10.set_xlabel(r'$p_Q$: strictness of isolation in $Q$', fontsize=14)
ax10.set_ylabel(r'$\sigma_{Q}$: relative frequency of testing in $Q$',\
                fontsize=14)
ax10.tick_params(axis="x", labelsize=14)
ax10.tick_params(axis="y", labelsize=14)
levels10 = [0.41, 0.42, 0.43, 0.44, 0.45, 0.46]
CS10 = ax10.contour(p_Q_heat, sigma_Q_heat, phistar_p_Q_sigma_Q,levels10,\
                    colors='w')
ax10.clabel(CS10, inline=1, fontsize=12)
cbar10 = fig2.colorbar(pl10, ax=ax10)
cbar10.ax.tick_params(labelsize=12)
ax10.text(-0.1, 1.15, string.ascii_uppercase[0], transform=ax10.transAxes, 
        size=20, weight='bold')

ax8.plot(sigma_U[::2],phistar_sigma_U[::2],'-',linewidth=2.5)
ax8.plot(sigma_U[::2],phistar_sigma_U_low_compl[::2],'--',linewidth=2.5)
ax8.set_xlabel(r'$\sigma_{U_2}$: relative frequency of testing in $U_2$',\
               fontsize=14)
ax8.set_xlim(np.min(sigma_U),np.max(sigma_U))
ax8.set_ylim(0.3,0.6)
ax8.tick_params(axis="x", labelsize=14)
ax8.tick_params(axis="y", labelsize=14)
ax8.legend(['strict quaran./isol. \n $(p_Q,p_I)=(0.2,0.1)$',\
            'weak quaran./isol. \n $(p_Q,p_I)=(0.7,0.35)$'],\
            prop=dict(size=12,weight='bold'))
ax8.text(-0.1, 1.15, string.ascii_uppercase[1], transform=ax8.transAxes, 
        size=20, weight='bold')

pl11 = ax11.pcolormesh(cov_heat, sigma_U_heat, phistar_cov_sigma_U, cmap=cmap,\
                       shading='auto')
ax11.axis([l_cov_heat, r_cov_heat, l_sigma_U_heat, r_sigma_U_heat])
ax11.set_ylim(np.min(sigma_U_heat),np.max(sigma_U_heat))
ax11.set_xlabel(r'$\omega$: tracing coverage', fontsize=14)
ax11.set_ylabel(r'$\sigma_{U_2}$: relative frequency of testing in $U_2$',\
                fontsize=14)
levels11 = [0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
CS11 = ax11.contour(cov_heat, sigma_U_heat, phistar_cov_sigma_U,levels11,\
                    colors='w')
ax11.clabel(CS11, inline=1, fontsize=12) 
ax11.tick_params(axis="x", labelsize=14)
ax11.tick_params(axis="y", labelsize=14)
cbar11 = fig2.colorbar(pl11, ax=ax11)
cbar11.ax.tick_params(labelsize=12)
ax11.text(-0.1, 1.15, string.ascii_uppercase[2], transform=ax11.transAxes, 
          size=20, weight='bold') 

#plt.savefig('Figure7.eps', format='eps', dpi=600)
#plt.savefig('Figure7.pdf', format='pdf', dpi=600)




    
  
    
  
