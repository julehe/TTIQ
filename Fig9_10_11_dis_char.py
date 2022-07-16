# Figure 9 & 10 & 11: phistar varying disease characteristics

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

# baseline parameter setting 
cov_b = 0.65 # tracing coverage
sigma_U_b = 93 # rel. freq. of testing in U2
sigma_Q_b = 300 # rel. freq. of testing in Q
p_I_b = 0.1 # strictness of isolation
p_Q_b = 0.2 # strictness of quarantine
tau_b = 2 # tracing delay
theta_b = 1.5 # scaling factor of presymptomatic transmission rate
lat_b = 3.5 # average duraiton of latency period
pre_b = 2 # average duration of presymptomatic infectious period
inf_b = 9 # average total duration of infectious period
R0 = 3.3 # basic reproduction number

# mesh-finess and colormap for heatmaps
n = 100
cmap = plt.get_cmap('viridis')

###################
### theta & pre ###
###################

theta, pre = np.meshgrid(np.linspace(0.5,2.5,n), np.linspace(1,4,n))
phistar_theta_pre = []

for o in range(0,len(theta)):
    entry = [];
    thetas = theta[1];
    pres = pre[:,1];
    for j in range(0,len(pre)):
        entry.append(phistar([cov_b, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b,\
                              tau_b, thetas[j], lat_b, pres[o], inf_b, R0],0))
    phistar_theta_pre.append(entry)
    
phistar_theta_pre = np.array(phistar_theta_pre)
l_theta=theta.min()
r_theta=theta.max()
l_pre=pre.min()
r_pre=pre.max()
l_c_theta_pre,r_c_theta_pre  = np.abs(phistar_theta_pre).min(),\
                               np.abs(phistar_theta_pre).max()

#################
### lat & inf ###
#################

inf, lat = np.meshgrid(np.linspace(6,12,n), np.linspace(1,6,n))
phistar_lat_inf = []

for o in range(0,len(inf)):
    entry = [];
    infs = inf[1];
    lats = lat[:,1];
    for j in range(0,len(lat)):
        entry.append(phistar([cov_b, sigma_U_b, sigma_Q_b, p_I_b, p_Q_b,\
                              tau_b, theta_b, lats[o], pre_b, infs[j], R0],0))
    phistar_lat_inf.append(entry)
    
phistar_lat_inf = np.array(phistar_lat_inf)
l_inf=inf.min()
r_inf=inf.max()
l_lat=lat.min()
r_lat=lat.max()
l_c_lat_inf,r_c_lat_inf  = np.abs(phistar_lat_inf).min(), np.abs(phistar_lat_inf).max()

##########################################################################
### phistar without any TTIQ (no TTIQ) vs. only testing (only testing) ### 
### vs. testing and tracing (full TTIQ) for different lat & inf        ###
##########################################################################

# lat = 2, inf =7
no_TTI_2_7 = phistar([cov_b,sigma_U_b,sigma_Q_b,1,1,tau_b,theta_b,2,pre_b,\
                      7,R0],0)

only_testing_2_7 = phistar([0,sigma_U_b,sigma_Q_b,0.1,0.2,tau_b,theta_b,2,pre_b,\
                            7,R0],0)
full_TTI_2_7 = phistar([cov_b,sigma_U_b,sigma_Q_b,0.1,0.2,tau_b,theta_b,2,pre_b,\
                        7,R0],0)
phistar_2_7 = [no_TTI_2_7, only_testing_2_7, full_TTI_2_7]

# lat = 5, inf = 7
no_TTI_5_7 = phistar([cov_b,sigma_U_b,sigma_Q_b,1,1,tau_b,theta_b,5,pre_b,\
                      7,R0],0)
only_testing_5_7 = phistar([0,sigma_U_b,sigma_Q_b,0.1,0.2,tau_b,theta_b,5,pre_b,\
                            7,R0],0)
full_TTI_5_7 = phistar([cov_b,sigma_U_b,sigma_Q_b,0.1,0.2,tau_b,theta_b,5,pre_b,\
                        7,R0],0)
phistar_5_7 = [no_TTI_5_7, only_testing_5_7, full_TTI_5_7]

# lat = 2, inf = 11
no_TTI_2_11 = phistar([cov_b,sigma_U_b,sigma_Q_b,1,1,tau_b,theta_b,2,pre_b,\
                       11,R0],0)
only_testing_2_11 = phistar([0,sigma_U_b,sigma_Q_b,0.1,0.2,tau_b,theta_b,2,pre_b,\
                             11,R0],0)
full_TTI_2_11 = phistar([cov_b,sigma_U_b,sigma_Q_b,0.1,0.2,tau_b,theta_b,2,pre_b,\
                         11,R0],0)
phistar_2_11 = [no_TTI_2_11, only_testing_2_11, full_TTI_2_11]
    
# lat = 5, inf = 11
no_TTI_5_11 = phistar([cov_b,sigma_U_b,sigma_Q_b,1,1,tau_b,theta_b,5,pre_b,\
                       11,R0],0)
only_testing_5_11 = phistar([0,sigma_U_b,sigma_Q_b,0.1,0.2,tau_b,theta_b,5,pre_b,\
                             11,R0],0)
full_TTI_5_11 = phistar([cov_b,sigma_U_b,sigma_Q_b,0.1,0.2,tau_b,theta_b,5,pre_b,\
                         11,R0],0)
phistar_5_11 = [no_TTI_5_11, only_testing_5_11, full_TTI_5_11]

names = ["no TTIQ", "only testing", "full TTIQ"]

#####################################
### (cov & tau) for different lat ###
#####################################

# baseline duration of latent phase
tau, cov = np.meshgrid(np.linspace(0.5, 14, n), np.linspace(0, 1, n))
phistar_cov_tau_lat_bas = []
for o in range(0,n):
    entry = [];
    taus = tau[1];
    covs = cov[:,1];
    for j in range(0,n):
        entry.append(phistar([covs[o], sigma_U_b, sigma_Q_b, p_I_b, p_Q_b,\
                              taus[j], theta_b, lat_b, pre_b, inf_b, R0],0))
    phistar_cov_tau_lat_bas.append(entry)
    
phistar_cov_tau_lat_bas = np.array(phistar_cov_tau_lat_bas)
l_tau=tau.min()
r_tau=tau.max()
l_cov=cov.min()
r_cov=cov.max()
l_c_cov_del1,r_c_cov_del1  = np.abs(phistar_cov_tau_lat_bas).min(),\
                             np.abs(phistar_cov_tau_lat_bas).max()

# longer duration of latent phase
phistar_cov_tau_lat_long = []
for o in range(0,n):
    entry = [];
    taus = tau[1];
    covs = cov[:,1];
    for j in range(0,n):
        entry.append(phistar([covs[o], sigma_U_b, sigma_Q_b, p_I_b, p_Q_b,\
                              taus[j], theta_b, 10, pre_b, inf_b, R0],0))
    phistar_cov_tau_lat_long.append(entry)
    
phistar_cov_tau_lat_long = np.array(phistar_cov_tau_lat_long)
l_c_cov_del2,r_c_cov_del2  = np.abs(phistar_cov_tau_lat_long).min(),\
                             np.abs(phistar_cov_tau_lat_long).max()

##############
### FIGURE ###
##############

fig1 = plt.figure(figsize=(5.5, 4),constrained_layout=True)
grid = fig1.add_gridspec(ncols=1, nrows=1, wspace=0.1)
ax1 = fig1.add_subplot(grid[0])

pl1 = ax1.pcolormesh(theta, pre, phistar_theta_pre, cmap=cmap,\
                     vmin=l_c_theta_pre, vmax=r_c_theta_pre, shading='auto')
ax1.axis([l_theta, r_theta, l_pre, r_pre])
ax1.set_xlabel(r'$\theta$: scaling factor presym. transmissions', fontsize=14)
ax1.set_ylabel(r'$1/\gamma_1$: av. duration presym. phase', fontsize=14)
ax1.tick_params(axis="x", labelsize=14)
ax1.tick_params(axis="y", labelsize=14)
levels1 = [0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54]
CS1 = ax1.contour(theta, pre, phistar_theta_pre, levels1, colors='w')
ax1.clabel(CS1, inline=1, fontsize=12)
cbar1 = fig1.colorbar(pl1, ax=ax1)
cbar1.ax.tick_params(labelsize=14)

#plt.savefig('Figure9.pdf', format='pdf', dpi=600)
#plt.savefig('Figure9.eps', format='eps', dpi=600)

fig2 = plt.figure(figsize=(14, 8),constrained_layout=True)
grid = fig2.add_gridspec(ncols=3, nrows=2, wspace=0.1)
ax2 = fig2.add_subplot(grid[0])
ax3 = fig2.add_subplot(grid[1])
ax4 = fig2.add_subplot(grid[2])
ax5 = fig2.add_subplot(grid[4])
ax6 = fig2.add_subplot(grid[5])

pl2 = ax2.pcolormesh(inf, lat, phistar_lat_inf, cmap=cmap, vmin=l_c_lat_inf,\
                     vmax=r_c_lat_inf, shading='auto')
ax2.axis([l_inf, r_inf, l_lat, r_lat])
ax2.set_xlabel(r'$1/\gamma$: av. duration infectious phase', fontsize=14)
ax2.set_ylabel(r'$1/\alpha$: av. duration latent phase', fontsize=14)
ax2.tick_params(axis="x", labelsize=14)
ax2.tick_params(axis="y", labelsize=14)
levels2 = [0.4, 0.42, 0.44, 0.46, 0.48, 0.5]
CS2 = ax2.contour(inf, lat, phistar_lat_inf, levels2, colors='w')
ax2.clabel(CS2, inline=1, fontsize=12)
cbar2 = fig2.colorbar(pl2, ax=ax2)
cbar2.ax.tick_params(labelsize=14) 
ax2.text(-0.1, 1.15, string.ascii_uppercase[0], transform=ax2.transAxes,
        size=20, weight='bold')

ax3.bar(names,phistar_2_7 ,width=0.4,color=[My_Blue])
ax3.set_ylim(0,0.6)
ax3.tick_params(axis="x", labelsize=14)
ax3.tick_params(axis="y", labelsize=14)
ax3.set_title(r'$1/\alpha=2, 1/\gamma=7$', fontsize=16)
ax3.annotate(str(round(phistar_2_7 [0],3)), xy=(names[0],phistar_2_7 [0]),\
             ha='center', va='bottom', fontsize=14)
ax3.annotate(str(round(phistar_2_7 [1],3)), xy=(names[1],phistar_2_7 [1]),\
             ha='center', va='bottom', fontsize=14)
ax3.annotate(str(round(phistar_2_7 [2],3)), xy=(names[2],phistar_2_7 [2]),\
             ha='center', va='bottom', fontsize=14)
ax3.text(-0.1, 1.15, string.ascii_uppercase[1], transform=ax3.transAxes,
        size=20, weight='bold')

ax4.bar(names,phistar_5_7,width=0.4,color=[My_Blue])
ax4.set_ylim(0,0.6)
ax4.tick_params(axis="x", labelsize=14)
ax4.tick_params(axis="y", labelsize=14)
ax4.set_title(r'$1/\alpha=5, 1/\gamma=7$', fontsize=16)
ax4.annotate(str(round(phistar_5_7[0],3)), xy=(names[0],phistar_5_7[0]),\
             ha='center', va='bottom', fontsize=14)
ax4.annotate(str(round(phistar_5_7[1],3)), xy=(names[1],phistar_5_7[1]),\
             ha='center', va='bottom', fontsize=14)
ax4.annotate(str(round(phistar_5_7[2],3)), xy=(names[2],phistar_5_7[2]),\
             ha='center', va='bottom', fontsize=14)
ax4.text(-0.1, 1.15, string.ascii_uppercase[2], transform=ax4.transAxes,
        size=20, weight='bold')

ax5.bar(names,phistar_2_11,width=0.4,color=[My_Blue])
ax5.set_ylim(0,0.6)
ax5.tick_params(axis="x", labelsize=14)
ax5.tick_params(axis="y", labelsize=14)
ax5.set_title(r'$1/\alpha=2, 1/\gamma=11$', fontsize=16)
ax5.annotate(str(round(phistar_2_11[0],3)), xy=(names[0],phistar_2_11[0]),\
             ha='center', va='bottom', fontsize=14)
ax5.annotate(str(round(phistar_2_11[1],3)), xy=(names[1],phistar_2_11[1]),\
             ha='center', va='bottom', fontsize=14)
ax5.annotate(str(round(phistar_2_11[2],3)), xy=(names[2],phistar_2_11[2]),\
             ha='center', va='bottom', fontsize=14)
ax5.text(-0.1, 1.15, string.ascii_uppercase[3], transform=ax5.transAxes,
        size=20, weight='bold')

ax6.bar(names,phistar_5_11,width=0.4,color=[My_Blue])
ax6.set_ylim(0,0.6)
ax6.tick_params(axis="x", labelsize=14)
ax6.tick_params(axis="y", labelsize=14)
ax6.set_title(r'$1/\alpha=5, 1/\gamma=11$', fontsize=16)
ax6.annotate(str(round(phistar_5_11[0],3)), xy=(names[0],phistar_5_11[0]),\
             ha='center', va='bottom', fontsize=14)
ax6.annotate(str(round(phistar_5_11[1],3)), xy=(names[1],phistar_5_11[1]),\
             ha='center', va='bottom', fontsize=14)
ax6.annotate(str(round(phistar_5_11[2],3)), xy=(names[2],phistar_5_11[2]),\
             ha='center', va='bottom', fontsize=14)
ax6.text(-0.1, 1.15, string.ascii_uppercase[4], transform=ax6.transAxes,
        size=20, weight='bold')

#plt.savefig('Figure10.pdf', format='pdf', dpi=600)
#plt.savefig('Figure10.eps', format='eps', dpi=600)

fig3 = plt.figure(figsize=(9.9, 3.96),constrained_layout=True)
grid = fig3.add_gridspec(ncols=2, nrows=1, wspace=0.1)
ax8 = fig3.add_subplot(grid[0])
ax9 = fig3.add_subplot(grid[1])

pl8 = ax8.pcolormesh(tau, cov, phistar_cov_tau_lat_bas, cmap=cmap,\
                     vmin=l_c_cov_del1, vmax=r_c_cov_del1, shading='auto')
ax8.axis([l_tau, r_tau, l_cov, r_cov])
ax8.set_xlabel(r'$\kappa$: tracing delay', fontsize=14)
ax8.set_ylabel(r'$\omega$: tracing coverage', fontsize=14)
ax8.set_title('average latent phase $3.5$ days', fontsize=16)
ax8.tick_params(axis="x", labelsize=14)
ax8.tick_params(axis="y", labelsize=14)
levels8 = [0.41, 0.42, 0.44, 0.46, 0.48, 0.5]
CS8 = ax8.contour(tau, cov, phistar_cov_tau_lat_bas, levels8, colors='w')
ax8.clabel(CS8, inline=1, fontsize=12)
cbar8 = fig3.colorbar(pl8, ax=ax8)
cbar8.ax.tick_params(labelsize=14) 
ax8.text(-0.1, 1.15, string.ascii_uppercase[0], transform=ax8.transAxes,
        size=20, weight='bold')

pl9 = ax9.pcolormesh(tau, cov, phistar_cov_tau_lat_long, cmap=cmap,\
                     vmin=l_c_cov_del2, vmax=r_c_cov_del2, shading='auto')
ax9.axis([l_tau, r_tau, l_cov, r_cov])
ax9.set_xlabel(r'$\kappa$: tracing delay', fontsize=14)
ax9.set_ylabel(r'$\omega$: tracing coverage', fontsize=14)
ax9.set_title('average latent phase $10$ days', fontsize=16)
ax9.tick_params(axis="x", labelsize=14)
ax9.tick_params(axis="y", labelsize=14)
levels9 = [0.41, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54]
CS9 = ax9.contour(tau, cov, phistar_cov_tau_lat_long, levels9, colors='w')
ax9.clabel(CS9, inline=1, fontsize=12)
cbar9 = fig3.colorbar(pl9, ax=ax9, ticks=[0.42, 0.44, 0.46, 0.48, 0.50, 0.52,\
                                          0.54])
cbar9.ax.tick_params(labelsize=14) 
ax9.text(-0.1, 1.15, string.ascii_uppercase[1], transform=ax9.transAxes,
        size=20, weight='bold')

#plt.savefig('Figure11.pdf', format='pdf', dpi=600)
#plt.savefig('Figure11.eps', format='eps', dpi=600)
