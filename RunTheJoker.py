#imports
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


DATA_PATH = os.getenv("DATA_PATH", "/users/EllaMathews/Research/") #environment variable 
OUTPUT_PATH = os.getenv("OUTPUT_PATH", "/users/EllaMathews/Summer-Research") 

#random generator to ensure reproducibility? - from the joker tutorial
rnd = np.random.default_rng(seed=42)

#importing new data
new_6866 = QTable.read(f'{DATA_PATH}/rcat_ngc6866_v0.fits')
new_6811 = QTable.read(f'{DATA_PATH}/rcat_ngc6811_v0.fits')

def RunTheJoker(id_num, num_priors, mpi = False):
    new_ids_6811 = new_6811['GAIAEDR3_ID']
    new_ids_6866 = new_6866['GAIAEDR3_ID']

    datamatched6811 = new_6811[id_num == new_ids_6811]
    datamatched6866 = new_6866[id_num == new_ids_6866]

    matched = vstack([datamatched6811, datamatched6866])
    
    if len(matched) == 0:
        print("No RV data for this ID")
        return
        
    t1 = Time(matched["DATE-OBS"], format = "fits", scale = "tcb")
    data = tj.RVData(t = t1, rv = matched['vrad']*(u.kilometer/u.second), rv_err = matched['vrad_err']*(u.kilometer/u.second)) 

    if len(data) < 3:
        print("Not enough RV data")
        return

    if os.path.exists(f'{OUTPUT_PATH}/{id_num}') == False:
        os.makedirs(f'{OUTPUT_PATH}/{id_num}')
        os.makedirs(f'{OUTPUT_PATH}/{id_num}/Plots')

    fig1, ax1 = plt.subplots()
    _ = data.plot() #plotting rv vs time
    fig1.savefig(f"{OUTPUT_PATH}/{id_num}/Plots/RVvsTime_{id_num}") #saving figure to plots folder in research folder
    
    prior = tj.JokerPrior.default( #initializing the default prior
        P_min = 2 * u.day,
        P_max = 1e3 * u.day,
        sigma_K0 = 30 * u.km / u.s,
        sigma_v = 100 * u.km / u.s,
    )

    if os.path.exists(f"{OUTPUT_PATH}/{id_num}/prior_samples_{id_num}.hdf5"):
        prior_samples = tj.JokerSamples.read(f"{OUTPUT_PATH}/{id_num}/prior_samples_{id_num}.hdf5")
    else:
        prior_samples = prior.sample(size = num_priors, rng = rnd) #generating prior samples
        prior_samples.write(f"{OUTPUT_PATH}/{id_num}/prior_samples_{id_num}.hdf5", overwrite = True) #write out prior samples to research folder 

    if mpi is True:
        with schwimmbad.MultiPool() as pool:
            print("Multiprocessing")
            try:
                joker = tj.TheJoker(prior, rng=rnd, pool=pool)
                joker_samples = joker.rejection_sample(data, prior_samples, max_posterior_samples=256)
                print("done sampling")
            except:
                print("failed")
                return
    else:
        pool = None
        joker = tj.TheJoker(prior, rng=rnd) #creating instance of The Joker
        joker_samples = joker.rejection_sample(data, prior_samples, max_posterior_samples=256) #creating rejection samples 
   
    joker_samples.write(f"{OUTPUT_PATH}/{id_num}/rejection_samples_{id_num}.hdf5", overwrite = True) #writing out posterior samples (not MCMC)

    fig2, ax2 = plt.subplots()
    _ = tj.plot_rv_curves(joker_samples, data=data) #plotting RV curves from rejection sampler
    fig2.savefig(f"{OUTPUT_PATH}/{id_num}/Plots/RVCurves_{id_num}") #saving figure to plots folder in research folder

    #plotting period against eccentricity
    fig3, ax3 = plt.subplots()
    with quantity_support():
        ax3.scatter(joker_samples["P"], joker_samples["e"], s=20, lw=0, alpha=0.5)
    ax3.set_xscale("log")
    ax3.set_xlim(1, 1e3)
    ax3.set_ylim(0, 1)
    ax3.set_xlabel("$P$ [day]")
    ax3.set_ylabel("$e$")
    fig3.savefig(f"{OUTPUT_PATH}/{id_num}/Plots/PeriodvsEccent_{id_num}") #saving figure to plots folder in research folder 

    if len(joker_samples) == 1: #if only one sample need MCMC
        #MCMC with NUTS sampler 
        with prior.model:
            mcmc_init = joker.setup_mcmc(data, joker_samples)
            trace = pm.sample(tune=500, draws=500, start=mcmc_init, cores=1, chains=2)
        mcmc_samples = tj.JokerSamples.from_inference_data(prior, trace, data) #convert trace into jokersamples
        mcmc_samples.write(f'{OUTPUT_PATH}/{id_num}/rejection_samples_MCMC_{id_num}.hdf5', overwrite = True) #write out MCMC posterior samples 
        
        fig4, ax4 = plt.subplots()
        _ = tj.plot_rv_curves(mcmc_samples, data=data) #plotting RV curves from MCMC rejection sampler
        fig4.savefig(f"{OUTPUT_PATH}/{id_num}/Plots/RVCurves_MCMC_{id_num}") #saving figure to plots folder in research folder
        
        #plotting period vs eccentricity
        fig5, ax5 = plt.subplots()
        with quantity_support():
            ax5.scatter(mcmc_samples["P"], mcmc_samples["e"], s=20, lw=0, alpha=0.5)
        ax5.set_xscale("log")
        ax5.set_xlim(1, 1e3)
        ax5.set_ylim(0, 1)
        ax5.set_xlabel("$P$ [day]")
        ax5.set_ylabel("$e$")
        fig5.savefig(f"{OUTPUT_PATH}/{id_num}/Plots/PeriodvsEccent_MCMC_{id_num}") #saving figure to plots folder in research folder
        
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('id', help = 'star id', type = int)
    parser.add_argument('--prior', help = 'num of prior samp default 1000', type = int, default = 1000)
    args = parser.parse_args()
    
    RunTheJoker(args.id, args.prior)
