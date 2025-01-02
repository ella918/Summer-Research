from __future__ import division, print_function, absolute_import

#matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.ascii as at
import astropy.units as u
import astropy.table as astropy
import astropy.units as u
import matplotlib.pyplot as plt 
import numpy as np 
import thejoker as tj
import h5py 
from astropy.table import QTable, Table, Column 
from astropy.time import Time
from astropy.visualization.units import quantity_support 
import os 


workpath2 = '/data2/labs/douglste-laf-lab/mathewea/50.0M'
unimodal_table2 = Table.read('/data2/labs/douglste-laf-lab/mathewea/Summer-Research/unimodalcheck_rerun.csv')


def circ_function(Porb,Pcirc,alpha=0.35,beta=0.14,gamma=1.0):
    """
    Compute the eccentricity distribution for the input
    Porb values. Alpha, beta, and gamma are taken from 
    Meibom & Mathieu (2005). 
    
    Inputs:
        Porb: array-like (Quantity optional)
        Pcirc: float (Quantity optional)
        alpha, beta, gamma: floats (optional)
        
    Outputs:
        eccentricity: array
    """
    
    eccentricities = np.zeros_like(Porb)
    
    # If Porb <= Pcirc, then e=0, so we don't need to change anything
    
    # If Porb > Pcirc, compute the circularization function
    
    gtr_pc = Porb > Pcirc
    
    cf_part = np.e**(beta * (Pcirc - Porb[gtr_pc]))
    
    eccentricities[gtr_pc] = alpha * (1 - cf_part)**gamma
    
    return eccentricities
ptest = np.logspace(0,3,100)
etest = circ_function(ptest,10.2)
# m35 = at.read("/home/stephanie/data/catalogs/M35_orbits_meibommathieu2005.csv")

# import emcee
# emcee.__version__

# def lnlike(Pcirc,Porb,obs_ecc,obs_err):
#     mod_ecc = circ_function(Porb,Pcirc)
    
#     sigma2 = obs_err**2 + mod_ecc**2
#     inv_sigma2 = 1 / sigma2
#     return -0.5*np.sum((obs_ecc - mod_ecc)**2 * inv_sigma2
#                       - np.log(inv_sigma2))
# def lnprior(Pcirc):
#     if Pcirc>=0.01 and Pcirc<50:
#         return 0.0
#     else:
#         return -np.inf

# def lnprob(Pcirc,Porb,obs_ecc,obs_err):
    
#     lp = lnprior(Pcirc)
#     if not np.isfinite(lp):
#         return - np.inf
#     else:
#         return lp + lnlike(Pcirc,Porb,obs_ecc,obs_err)

plt.figure()
ax = plt.subplot(111)
plt.plot(np.log10(ptest),etest)


data_for_plots = QTable()
ids = []
P_median = []
P_lower = []
P_upper = []
e_median = []
e_lower = []
e_upper = []
MCMC = []
#MAP = []

for i in range(len(unimodal_table2)):
    if unimodal_table2['unimodal'][i] == 1:
        idnum = unimodal_table2['id'][i]
        ids.append(idnum)
        if unimodal_table2['MCMC'][i] == 0:
            joker_samples = tj.JokerSamples.read(f'{workpath2}/{idnum}/rejection_samples_50.0M_{idnum}.hdf5')
            MCMC.append(0)

        if unimodal_table2['MCMC'][i] == 1:
            joker_samples = tj.JokerSamples.read(f'{workpath2}/{idnum}/rejection_samples_MCMC_50.0M_{idnum}.hdf5')
            MCMC.append(1)

        p_median = np.percentile(joker_samples['P'], 50)
        #print(p_median)
        p_16 = np.percentile(joker_samples['P'], 16)
        #print(p_16)
        p_84 = np.percentile(joker_samples['P'], 84)
        #print(p_84)

        e_median1 = np.percentile(joker_samples['e'], 50)
        #print(e_median1)
        e_16 = np.percentile(joker_samples['e'], 16)
        #print(e_16)
        e_84 = np.percentile(joker_samples['e'], 84)
        #print(e_84)

        P_median.append(p_median)
        P_lower.append(p_median - p_16)
        P_upper.append(p_84 - p_median)

        e_median.append(e_median1)
        e_lower.append(e_median1 - e_16)
        e_upper.append(e_84 - e_median1)

        #Map = tj.MAP_sample(joker_samples)
        #MAP.append(Map)

data_for_plots['id'] = ids
data_for_plots['MCMC'] = MCMC
data_for_plots['P_median'] = P_median
data_for_plots['P_lower'] = P_lower
data_for_plots['P_upper'] = P_upper 
data_for_plots['e_median'] = e_median
data_for_plots['e_lower'] = e_lower
data_for_plots['e_upper'] = e_upper 
# #data_for_plots['MAP'] = MAP 

data_for_plots.write('data_for_unimodal_plots_reran.csv', format = 'csv', overwrite = True)
print('made table')


#plotting e vs P plot with uncertainties 
x = np.array(data_for_plots['P_median'])
y = np.array(data_for_plots['e_median'])
asymmetric_error_x = np.array([np.array(data_for_plots['P_lower']), np.array(data_for_plots['P_upper'])])
#print(asymmetric_error_x)
asymmetric_error_y = np.array([np.array(data_for_plots['e_lower']), np.array(data_for_plots['e_upper'])])
#print(asymmetric_error_y)
ax.errorbar(x, y, xerr = asymmetric_error_x, yerr = asymmetric_error_y, fmt = 'o')
ax.set_xlabel('P (d)')
ax.set_xscale('log')
ax.set_ylabel('e')
ax.set_title('e vs P for Unimodal Stars in NGC6811 and NGC6866')
plt.show()
plt.savefig('Unimodal_Plots_AAS.png')
