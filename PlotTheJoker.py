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

def PlotTheJoker(id_num):

	if os.path.exists(f"{DATA_PATH}/{id_num}/prior_samples_{id_num}.hdf5") is False:
		print("Not enough RV Data")
		return

	prior_samples = tj.JokerSamples.read(f"{DATA_PATH}/{id_num}/prior_samples_{id_num}.hdf5")
	joker_samples = tj.JokerSamples.read(f"{id_num}/rejection_samples_{id_num}.hdf5")


	fig1, ax1 = plt.subplots()
    _ = tj.plot_rv_curves(joker_samples, data=data) #plotting RV curves from rejection sampler
    fig1.savefig(f"{id_num}/RVCurves_{id_num}, overwrite = True") #saving figure to plots folder in research folder
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
    fig2.savefig(f"{id_num}/PeriodvsEccent_{id_num}", overwrite = True) #saving figure to plots folder in research folder 
    print("Period vs Eccentricity plotted")

    if len(joker_samples) == 1: 
    	mcmc_samples = tj.JokerSamples.read(f'{id_num}/rejection_samples_MCMC_{id_num}.hdf5')

    	fig3, ax3 = plt.subplots()
        _ = tj.plot_rv_curves(mcmc_samples, data=data) #plotting RV curves from MCMC rejection sampler
        fig3.savefig(f"{id_num}/RVCurves_MCMC_{id_num}", overwrite = True) #saving figure to plots folder in research folder
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
        fig4.savefig(f"{id_num}/PeriodvsEccent_MCMC_{id_num}", overwrite = True) #saving figure to plots folder in research folder
        print("Period vs Eccentricity from MCMC plotted")

  	return
