"""
Figure 7: global sensitivity of $\phi^*$ to changes in TTIQ parameters using
latin hypercube samping (LHS) and partial rank correlation coefficients 
(PRCCs).
"""

import numpy as np
from phistar import phistar
import scipy as sc
from scipy.stats import qmc
import pingouin as pg
import pandas as pd
import matplotlib.pyplot as plt

# costum colors for plotting
My_Amber = "#E69F00"
My_Blue = "#56B4E9"
My_Green = "#009E73"
My_Violet = '#661D98'

color_list = [My_Blue, My_Amber, My_Green, My_Violet]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)

######################################
### LHS and calculation of phistar ###
######################################

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

# number of samples
runs = 150000

sampler = qmc.LatinHypercube(d=5)
sample = sampler.random(n=runs)

# lower and upper bounds for intervals to sample TTIQ parameter from
# cov, sigma_U, sigma_Q, p_I, tau
l_bounds = [0, 1, 1, 0, 0.5]
u_bounds = [1, 186, 600, 1, 14]

# scale samples to bounds
# cov, sigma_U, sigma_Q, p_I, tau
sample = qmc.scale(sample, l_bounds, u_bounds)

# calulate phistar corresponding to samples
outcome = []
for k in range(0,runs):
    outcome.append(phistar([sample[k,0],sample[k,1],sample[k,2],sample[k,3],3,\
                            sample[k,4],theta_b,lat_b,pre_b,inf_b,R0],1))

###################################
### rank transformation of data ###
##################################

for k in range(0,5):
    sample[:,k]=sc.stats.rankdata(sample[:,k], method='average')

outcome_ranked = sc.stats.rankdata(outcome, method='average')
outcome_ranked.shape = (runs,1)

all_data = np.append(sample,outcome_ranked,1)
all_data = pd.DataFrame(data=all_data, columns=["cov", "sigma_U", "sigma_Q",\
                                                "p_I","tau","phi"])
    
############################
### calculation of PRCCs ###
############################

pear_cov = pg.partial_corr(data=all_data, x='cov', y='phi',\
                           covar=["sigma_U","sigma_Q","p_I","tau"])
r_cov = pear_cov.loc['pearson','r']

pear_sigma_U = pg.partial_corr(data=all_data, x='sigma_U', y='phi',\
                               covar=["cov","sigma_Q","p_I","tau"])
r_sigma_U = pear_sigma_U.loc['pearson','r']

pear_sigma_Q = pg.partial_corr(data=all_data, x='sigma_Q', y='phi',\
                               covar=["cov","sigma_U","p_I","tau"])
r_sigma_Q = pear_sigma_Q.loc['pearson','r']

pear_p_I = pg.partial_corr(data=all_data, x='p_I', y='phi',\
                           covar=["cov","sigma_U", "sigma_Q","tau"])
r_p_I = pear_p_I.loc['pearson','r']

pear_tau = pg.partial_corr(data=all_data, x='tau', y='phi',\
                           covar=["cov","sigma_U", "sigma_Q","p_I"])
r_tau = pear_tau.loc['pearson','r']

rs = [r_cov, r_sigma_Q, r_tau, r_sigma_U, r_p_I]

xticks = [1,2,3,4,5]
params = [r"$\omega$", r"$\sigma_{Q}$",r"$\kappa$", r"$\sigma_{U_2}$", r"$p_I$"]

fig = plt.figure(figsize=(4.7, 4),constrained_layout=True)
grid = fig.add_gridspec(ncols=1, nrows=1, wspace=0.1)
ax21 = fig.add_subplot(grid[0])
ax21.bar(xticks,rs,color=My_Blue)
ax21.set_xticks(xticks)
ax21.set_xticklabels(params)
ax21.tick_params(axis="x", labelsize=12)
ax21.tick_params(axis="y", labelsize=12)
ax21.plot((0, 6), (0, 0), 'k-',linewidth=1)
ax21.annotate(str(round(rs[0],2)), xy=(xticks[0],rs[0]), ha='center',\
              va='bottom', fontsize=11)
ax21.annotate(str(round(rs[1],2)), xy=(xticks[1],rs[1]), ha='center',\
              va='bottom', fontsize=11)
ax21.annotate(str(round(rs[2],2)), xy=(xticks[2],rs[2]), ha='center',\
              va='top', fontsize=11)
ax21.annotate(str(round(rs[3],2)), xy=(xticks[3],rs[3]), ha='center',\
              va='bottom', fontsize=11)
ax21.annotate(str(round(rs[4],2)), xy=(xticks[4],rs[4]), ha='center',\
              va='top', fontsize=11)
ax21.set_xlim(0,6)
ax21.set_ylim(-1,1)
ax21.set_xlabel('TTIQ parameters', fontsize=14)
ax21.set_ylabel('PRCC for $\\phi^*$', fontsize=14)  

#plt.savefig('Figure7.eps', format='eps', dpi=600)
#plt.savefig('Figure7.pdf', format='pdf', dpi=600)

# print("PRCC: tracing coverage")
# print(pear_cov)
# print("PRCC: rel. freq. of test. in Q")
# print(pear_sigma_Q)
# print("PRCC: tracing delay")
# print(pear_tau)
# print("PRCC: strictness of isolation")
# print(pear_p_I)
# print("PRCC: rel. freq. of test. in U2")
# print(pear_sigma_U)

