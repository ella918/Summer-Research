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


DATA_PATH = os.getenv("DATA_PATH", "/users/EllaMathews/Summer-Research/") #environment variable 
jobid = os.getenv("SLURM_JOB_ID", "-9999")
if jobid != "-9990":
    workpath = "/scratch/ella_oneprior/"
else:
    workpath = DATA_PATH

#random generator to ensure reproducibility
rnd = np.random.default_rng(seed=42)

def CreatePriors(num_priors):
	mils = num_priors/1000000
	prior = tj.JokerPrior.default( #initializing the default prior
        P_min = 2 * u.day,
        P_max = 1e3 * u.day,
        sigma_K0 = 30 * u.km / u.s,
        sigma_v = 100 * u.km / u.s,
    )
	print('initialized default prior')

	print('Creating new priors')
	prior_samples = prior.sample(size = num_priors, rng = rnd, return_logprobs = True) #generating prior samples
	print("New priors created")
	prior_samples.write(f"{workpath}prior_samples_{mils:.0f}M.hdf5", overwrite = True) #write out prior samples to research folder 
	print("Priors saved to file!")
	return

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--prior', help = 'num of prior samp default 300000000', type = int, default = 300000000)
	args = parser.parse_args()
	print('args')

	CreatePriors(args.prior)

 









