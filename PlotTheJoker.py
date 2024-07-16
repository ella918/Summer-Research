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

def PlotTheJoker(id_num):

        if os.path.exists(f"{DATA_PATH}/{id_num}/prior_samples_{id_num}.hdf5") == False: #checking if the joker has been run 
                print("The Joker has not been run on this object yet.")
                return
        
        if os.path.exists(f"{DATA_PATH}/{id_num}/Plots") == False:
                os.makedirs(f"{DATA_PATH}/{id_num}/Plots") #creating a plots folder 
        #recreating the joker object to plot rv curves
        new_ids_6811 = new_6811['GAIAEDR3_ID']
        new_ids_6866 = new_6866['GAIAEDR3_ID']

        datamatched6811 = new_6811[id_num == new_ids_6811]
        datamatched6866 = new_6866[id_num == new_ids_6866]

        matched = vstack([datamatched6811, datamatched6866])
        t1 = Time(matched["DATE-OBS"], format = "fits", scale = "tcb")
        data = tj.RVData(t = t1, rv = matched['vrad']*(u.kilometer/u.second), rv_err = matched['vrad_err']*(u.kilometer/u.second))

        if len(matched) < 3:
                print("LESS THAN 3 DATA POINTS SOMETHING IS WRONG")

        #importing the outputs from running the joker (priors and  rejection samples)
        prior_samples = tj.JokerSamples.read(f"{DATA_PATH}/{id_num}/prior_samples_{id_num}.hdf5")
	joker_samples = tj.JokerSamples.read(f"{DATA_PATH}/{id_num}/rejection_samples_{id_num}.hdf5")

	#getting the lnk value to put on title (having trouble with this)
    K = joker_samples['K']
    K1st = np.percentile(K, 1)

	fig1, ax1 = plt.subplots()
	_ = tj.plot_rv_curves(joker_samples, data=data) #plotting RV curves from rejection sampler
	plt.title(f"KIC {id_num}, K1%={K1st}")
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
	plt.title(f"KIC {id_num}, K1%={K1st}") 
	print("Period vs Eccentricity plotted")

	if len(joker_samples) == 1: 
		mcmc_samples = tj.JokerSamples.read(f'{DATA_PATH}/{id_num}/rejection_samples_MCMC_{id_num}.hdf5')

		fig3, ax3 = plt.subplots()
		_ = tj.plot_rv_curves(mcmc_samples, data=data) #plotting RV curves from MCMC rejection sampler
		plt.title(f"KIC {id_num}, K1%={K1st}")
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
		plt.title(f"KIC {id_num}, K1%={K1st}")
	
		print("Period vs Eccentricity from MCMC plotted")

	return

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('id', help = 'star id', type = int)
	args = parser.parse_args()

	PlotTheJoker(args.id)
