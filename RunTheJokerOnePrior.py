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

DATA_PATH = os.getenv("DATA_PATH", "/users/EllaMathews/Summer-Research/") #environment variable 
jobid = os.getenv("SLURM_JOB_ID", "-9999")
if jobid != "-9990":
    workpath = "/scratch/ella_rerun/"
else:
    workpath = DATA_PATH
#random generator to ensure reproducibility
rnd = np.random.default_rng(seed=42)

#importing new data
new_6866 = QTable.read(f'{workpath}/rcat_ngc6866_v0.fits')
new_6811 = QTable.read(f'{workpath}/rcat_ngc6811_v0.fits')

def RunTheJokerOnePrior(id_num, mpi, num_priors):
    print(id_num in new_6866)
    print(id_num in new_6811)
    new_ids_6811 = new_6811['GAIAEDR3_ID']
    print(len(new_ids_6811))
    new_ids_6866 = new_6866['GAIAEDR3_ID']
    print(len(new_ids_6866))
    datamatched6811 = new_6811[id_num == new_ids_6811]
    print(len(datamatched6811))
    datamatched6866 = new_6866[id_num == new_ids_6866]
    print(len(datamatched6866))
    matched = vstack([datamatched6811, datamatched6866])
    print(len(matched))
    print("Created matched")
    if len(matched) == 0:
        print("No RV data for this ID")
        return

    os.makedirs(f'{workpath}/{id_num}')

    t1 = Time(matched["DATE-OBS"], format = "fits", scale = "tcb")
    data = tj.RVData(t = t1, rv = matched['vrad']*(u.kilometer/u.second), rv_err = matched['vrad_err']*(u.kilometer/u.second)) 
    print('created RV data object')

    mils = num_priors/1000000
    prior = tj.JokerPrior.default( #initializing the default prior
        P_min = 2 * u.day,
        P_max = 1e3 * u.day,
        sigma_K0 = 30 * u.km / u.s,
        sigma_v = 100 * u.km / u.s,
    )
    prior_samples = tj.JokerSamples.read(f'{workpath}/prior_samples_{mils}M.hdf5')

    if mpi is True: #multiprocessing
        with schwimmbad.MultiPool() as pool:
            print("Multiprocessing")
            try:
                joker = tj.TheJoker(prior, rng=rnd, pool=pool)
                joker_samples = joker.rejection_sample(data, prior_samples, max_posterior_samples=256, return_logprobs=True)
                print("done sampling")
            except:
                print("failed")
                #raise
                return
    else:
        pool = None
        joker = tj.TheJoker(prior, rng=rnd) #creating instance of The Joker
        joker_samples = joker.rejection_sample(data, prior_samples, max_posterior_samples=256, return_logprobs=True) #creating rejection samples 
    print('rejection sample created')
    joker_samples.write(f"{workpath}/{id_num}/rejection_samples_{mils}M.hdf5", overwrite = True) #writing out posterior samples (not MCMC)
    print('joker samples written out')


    if len(joker_samples) == 1: 
        print("1 sample, needs MCMC")#if only one sample need MCMC
        #MCMC with NUTS sampler 
        with prior.model:
            mcmc_init = joker.setup_mcmc(data, joker_samples)
            trace = pm.sample(tune=500, draws=500, start=mcmc_init, chains=2)
        mcmc_samples = tj.JokerSamples.from_inference_data(prior, trace, data) #convert trace into jokersamples
        mcmc_samples.write(f'{workpath}/{id_num}/rejection_samples_MCMC_{mils}M.hdf5', overwrite = True) #write out MCMC posterior samples 
    return 


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('id', help = 'star id', type = int)
    parser.add_argument('--mpi', help='False for no multiprocessing', type = bool, default = True)
    parser.add_argument('--prior', help = 'num of prior samp default 50000000', type = int, default = 50000000)
    args = parser.parse_args()
    print('args')

    RunTheJokerOnePrior(args.id, args.mpi, args.prior)