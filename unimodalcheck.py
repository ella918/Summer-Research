import astropy.table as at
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import h5py
import thejoker as tj
from astropy.table import QTable, Table, Column, vstack, unique
from astropy.time import Time
from astropy.visualization.units import quantity_support
import os


workpath = '/data2/labs/douglste-laf-lab/mathewea/200.0M_new'
idlist = QTable.read('/data2/labs/douglste-laf-lab/mathewea/Summer-Research/GAIADR3_IDs.csv')

new_6866 = QTable.read('/data2/labs/douglste-laf-lab/mathewea/rcat_ngc6866_v0.fits')
new_6811 = QTable.read('/data2/labs/douglste-laf-lab/mathewea/rcat_ngc6811_v0.fits')

datatable = Table()
ids = []
deltaT = []
mcmc = []
num_samples = []
unimodal = []
numRVs = []
maxgap = []
phasecov = []

for idnum in idlist['GAIAEDR3_ID']:
	new_ids_6811 = new_6811['GAIAEDR3_ID']
	new_ids_6866 = new_6866['GAIAEDR3_ID']
	datamatched6811 = new_6811[idnum == new_ids_6811]
	datamatched6866 = new_6866[idnum == new_ids_6866]
	matched = vstack([datamatched6811, datamatched6866])
	RV = len(matched)
	t1 = Time(matched["DATE-OBS"], format = "fits", scale = "tcb")
	data = tj.RVData(t = t1, rv = matched['vrad']*(u.kilometer/u.second), rv_err = matched['vrad_err']*(u.kilometer/u.second))

	if os.path.exists(f'{workpath}/{idnum}/rejection_samples_200.0M_{idnum}_new.hdf5'):
		ids.append(idnum)
		joker_samples = tj.JokerSamples.read(f'{workpath}/{idnum}/rejection_samples_200.0M_{idnum}_new.hdf5')
		numsamples = len(joker_samples)
		dt = t1[1:] - t1[:-1]
		if os.path.exists(f'{workpath}/{idnum}/rejection_samples_MCMC_200.0M_{idnum}_new.hdf5'):
			joker_samples = tj.JokerSamples.read(f'{workpath}/{idnum}/rejection_samples_MCMC_200.0M_{idnum}_new.hdf5')
			numsamples = 0
			mcmc_check = 1
			maxphase = tj.max_phase_gap(joker_samples, data)
			phasecov.append(maxphase)
			cov = tj.phase_coverage(joker_samples, data, n_bins=10)
			phasecov.append(cov)

		else:
			mcmc_check = 0
		if tj.is_P_unimodal(joker_samples, data):
			uni = 1
		else:
			uni = 0
		num_samples.append(numsamples)
		numRVs.append(RV)
		mcmc.append(mcmc_check)
		unimodal.append(uni)
		deltaT.append(dt.min().value)

datatable['id'] = ids
datatable['dt'] = deltaT
datatable['num_RV'] = numRVs
datatable['num_samples'] = num_samples
datatable['MCMC'] = mcmc
datatable['unimodal'] = unimodal
datatable['max_phase_gap'] = maxgap
datatable['phase_coverage'] = phasecov


datatable.write('unimodalcheck_200M_phase.csv', format = 'csv', overwrite = True)
