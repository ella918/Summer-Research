import astropy.table as at
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import h5py
import thejoker as tj
from astropy.table import QTable, Table, Column, vstack, unique
from astropy.time import Time
from astropy.visualization.units import quantity_support
import astropy.coordinates as coord
import pymc as pm
import arviz as az
import argparse
import os
import schwimmbad

DATA_PATH = os.getenv("DATA_PATH") 
workpath = f'/{DATA_PATH}/Summer-Research/'

#random generator to ensure reproducibility
rnd = np.random.default_rng(seed=42)

#importing new data
dataRV = QTable.read(f'{workpath}/MWE_RVData.csv')
print('read in data')
print(dataRV["DATE-OBS"])

id_num = 2128124963389008384
mpi = True
num_priors = 50000

t1 = Time(dataRV["DATE-OBS"], format = "isot", scale = "tcb")
data = tj.RVData(t = t1, rv = dataRV['vrad']*(u.kilometer/u.second), rv_err = dataRV['vrad_err']*(u.kilometer/u.second)) 
print('created RV data object')

mils = num_priors/1000000
print('mils')
prior = tj.JokerPrior.default( #initializing the default prior
    P_min = 2 * u.day,
    P_max = 1e3 * u.day,
    sigma_K0 = 30 * u.km / u.s,
    sigma_v = 100 * u.km / u.s,
)
prior_samples = prior.sample(size = num_priors, rng = rnd, return_logprobs=True)
prior_samples.write(f'{DATA_PATH}/prior_samples_{mils}M_{id_num}_MWE.hdf5', overwrite = True)

if mpi is True: #multiprocessing
    with schwimmbad.MultiPool() as pool:
        print("Multiprocessing")
        try:
            joker = tj.TheJoker(prior, rng=rnd, pool=pool)
            joker_samples = joker.rejection_sample(data, prior_samples, max_posterior_samples=256, return_logprobs=True)
            print("done sampling")
        except:
            print("failed")
            raise
            return
else:
    pool = None
    joker = tj.TheJoker(prior, rng=rnd) #creating instance of The Joker
    joker_samples = joker.rejection_sample(data, prior_samples, max_posterior_samples=256, return_logprobs=True) #creating rejection samples 
print('rejection sample created')
joker_samples.write(f"{DATA_PATH}/rejection_samples_{mils}M_{id_num}_MWE.hdf5", overwrite = True) #writing out posterior samples (not MCMC)
print('joker samples written out')
print(len(joker_samples), 'samples')

if len(joker_samples) == 1: 
    print("1 sample, needs MCMC")#if only one sample need MCMC
    #MCMC with NUTS sampler 
    with prior.model:
        mcmc_init = joker.setup_mcmc(data, joker_samples)
        trace = pm.sample(tune=500, draws=500, start=mcmc_init, chains=4)
    mcmc_samples = tj.JokerSamples.from_inference_data(prior, trace, data) #convert trace into jokersamples
    az.summary(trace, var_names=prior.par_names)
    az.plot_trace(trace, var_names = prior.par_names)
    plt.savefig(f'{DATA_PATH}/traceplot.png')
    mcmc_samples.write(f'{DATA_PATH}/rejection_samples_MCMC_{mils}M_{id_num}_MWE.hdf5', overwrite = True) #write out MCMC posterior samples  


#PLOTTING THE DATA

rvplotnoline = data.plot()
plt.savefig(f"{workpath}/RVvTime_{id_num}_MWE")

#getting the lnk value to put on title
K = joker_samples['K']
K1st = np.percentile(K, 1)
fig1, ax1 = plt.subplots()

_ = tj.plot_rv_curves(joker_samples, data=data) #plotting RV curves from rejection sampler
plt.title(f"ID: {id_num}, K1%={K1st}")
fig1.savefig(f"{workpath}/RVCurves_{id_num}_MWE") #saving figure to plots folder in script output folder
print("RV curves plotted")

#plotting period against eccentricity
fig2, ax2 = plt.subplots()
with quantity_support():
    ax2.scatter(joker_samples["P"], joker_samples["e"], s=20, lw=0, alpha=0.5)
ax2.set_xscale("log")
ax2.set_xlim(1, 1e3)
ax2.set_ylim(0, 1)
ax2.set_xlabel("$P$ [day]")
ax2.set_ylabel("$e$")
plt.title(f"ID: {id_num}, K1%={K1st}")

fig2.savefig(f"{workpath}/PeriodvsEccent_{id_num}_MWE") #saving figure to plots folder in script output  folder 
print("Period vs Eccentricity plotted")

if len(joker_samples) == 1: 

    fig3, ax3 = plt.subplots()
    _ = tj.plot_rv_curves(mcmc_samples, data=data) #plotting RV curves from MCMC rejection sampler
    plt.title(f"ID: {id_num}, K1%={K1st}")
    fig3.savefig(f"{workpath}/RVCurves_MCMC_{id_num}_MWE") #saving figure to plots folder in script output  folder
    print("RV curves from MCMC plotted")

    #plotting period vs eccentricity
    fig4, ax4 = plt.subplots()
    with quantity_support():
        ax4.scatter(mcmc_samples["P"], mcmc_samples["e"], s=20, lw=0, alpha=0.5)
    ax4.set_xscale("log")
    ax4.set_xlim(1, 1e3)
    ax4.set_ylim(0, 1)
    ax4.set_xlabel("$P$ [day]")
    ax4.set_ylabel("$e$")
    plt.title(f"ID: {id_num}, K1%={K1st}")
    fig4.savefig(f"{workpath}/PeriodvsEccent_MCMC_{id_num}_MWE") #saving figure to plots folder in script output  folder
    print("Period vs Eccentricity from MCMC plotted")
plt.close('all')
