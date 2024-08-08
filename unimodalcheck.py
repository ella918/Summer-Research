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



workpath2 = '/data2/labs/douglste-laf-lab/mathewea/TheJoker_Outputs'
workpath = '/data2/labs/douglste-laf-lab/mathewea/50.0M'

idlist2 = QTable.read('/data2/labs/douglste-laf-lab/mathewea/Summer-Research/GAIADR3_IDs.csv')
idlist = QTable.read('/data2/labs/douglste-laf-lab/mathewea/Summer-Research/rerunids.csv')
new_6866 = QTable.read('/data2/labs/douglste-laf-lab/mathewea/rcat_ngc6866_v0.fits')
new_6811 = QTable.read('/data2/labs/douglste-laf-lab/mathewea/rcat_ngc6811_v0.fits')

datatable = Table()
ids = []
mcmc = []
num_samples = []
unimodal = []
#bimodal = []
weirdids = [2075872941029111552, 2128018890577388672, 2128124963389008384, 2128145544872361984]
for idnum in idlist['id']:
	new_ids_6811 = new_6811['GAIAEDR3_ID']
	new_ids_6866 = new_6866['GAIAEDR3_ID']
	datamatched6811 = new_6811[idnum == new_ids_6811]
	datamatched6866 = new_6866[idnum == new_ids_6866]
	matched = vstack([datamatched6811, datamatched6866])
	t1 = Time(matched["DATE-OBS"], format = "fits", scale = "tcb")
	data = tj.RVData(t = t1, rv = matched['vrad']*(u.kilometer/u.second), rv_err = matched['vrad_err']*(u.kilometer/u.second))

	if os.path.exists(f'{workpath}/{idnum}/rejection_samples_50.0M_{idnum}.hdf5'):
		ids.append(idnum)
		joker_samples = tj.JokerSamples.read(f'{workpath}/{idnum}/rejection_samples_50.0M_{idnum}.hdf5')
		numsamples = len(joker_samples)
		if os.path.exists(f'{workpath}/{idnum}/rejection_samples_MCMC_50.0M_{idnum}.hdf5'):
			joker_samples = tj.JokerSamples.read(f'{workpath}/{idnum}/rejection_samples_MCMC_50.0M_{idnum}.hdf5')
			numsamples = 0
			mcmc_check = 1
		else:
			mcmc_check = 0
		if tj.is_P_unimodal(joker_samples, data):
			uni = 1
		else:
			uni = 0
		#if tj.is_P_Kmodal(joker_samples, data, 2):
		#	bi = 1
		#else:
		#	bi = 0 
		num_samples.append(numsamples)
		mcmc.append(mcmc_check)
		unimodal.append(uni)
		#bimodal.append(bi)

for idnum in weirdids:
	new_ids_6811 = new_6811['GAIAEDR3_ID']
	new_ids_6866 = new_6866['GAIAEDR3_ID']
	datamatched6811 = new_6811[idnum == new_ids_6811]
	datamatched6866 = new_6866[idnum == new_ids_6866]
	matched = vstack([datamatched6811, datamatched6866])
	t1 = Time(matched["DATE-OBS"], format = "fits", scale = "tcb")
	data = tj.RVData(t = t1, rv = matched['vrad']*(u.kilometer/u.second), rv_err = matched['vrad_err']*(u.kilometer/u.second))

	if os.path.exists(f'{workpath}/{idnum}/rejection_samples_50.0M_{idnum}.hdf5'):
		ids.append(idnum)
		joker_samples = tj.JokerSamples.read(f'{workpath}/{idnum}/rejection_samples_50.0M_{idnum}.hdf5')
		numsamples = len(joker_samples)
		if os.path.exists(f'{workpath}/{idnum}/rejection_samples_MCMC_50.0M_{idnum}.hdf5'):
			joker_samples = tj.JokerSamples.read(f'{workpath}/{idnum}/rejection_samples_MCMC_50.0M_{idnum}.hdf5')
			numsamples = 0
			mcmc_check = 1
		else:
			mcmc_check = 0
		if tj.is_P_unimodal(joker_samples, data):
			uni = 1
		else:
			uni = 0
		#if tj.is_P_Kmodal(joker_samples, data, 2):
		#	bi = 1
		#else:
		#	bi = 0 
		num_samples.append(numsamples)
		mcmc.append(mcmc_check)
		unimodal.append(uni)
		#bimodal.append(bi)

for idnum in idlist2['GAIAEDR3_ID']:
	if idnum in ids:
		continue

	new_ids_6811 = new_6811['GAIAEDR3_ID']
	new_ids_6866 = new_6866['GAIAEDR3_ID']
	datamatched6811 = new_6811[idnum == new_ids_6811]
	datamatched6866 = new_6866[idnum == new_ids_6866]
	matched = vstack([datamatched6811, datamatched6866])
	t1 = Time(matched["DATE-OBS"], format = "fits", scale = "tcb")
	data = tj.RVData(t = t1, rv = matched['vrad']*(u.kilometer/u.second), rv_err = matched['vrad_err']*(u.kilometer/u.second))

	if os.path.exists(f'{workpath2}/{idnum}/rejection_samples_{idnum}.hdf5'):
		ids.append(idnum)
		joker_samples = tj.JokerSamples.read(f'{workpath2}/{idnum}/rejection_samples_{idnum}.hdf5')
		numsamples = len(joker_samples)
		if os.path.exists(f'{workpath2}/{idnum}/rejection_samples_MCMC_{idnum}.hdf5'):
			joker_samples = tj.JokerSamples.read(f'{workpath2}/{idnum}/rejection_samples_MCMC_{idnum}.hdf5')
			numsamples = 0
			mcmc_check = 1
		else:
			mcmc_check = 0
		if tj.is_P_unimodal(joker_samples, data):
			uni = 1
		else:
			uni = 0
		#if tj.is_P_Kmodal(joker_samples, data, 2):
		#	bi = 1
		#else:
		#	bi = 0 
		num_samples.append(numsamples)
		mcmc.append(mcmc_check)
		unimodal.append(uni)
		#bimodal.append(bi)

datatable['id'] = ids
datatable['num_samples'] = num_samples
datatable['MCMC'] = mcmc
datatable['unimodal'] = unimodal
#datatable['bimodal'] = bimodal

ReRunTable = Table()
rerunids = []

datatable.write('unimodalcheck_rerun.csv', format = 'csv', overwrite = True)

# for i in range(len(datatable)):
# 	if datatable['MCMC'][i] == 0 and datatable['num_samples'][i] < 256:
# 		rerunids.append(datatable['id'][i])
# 	if datatable['MCMC'][i] == 1 and datatable['unimodal'][i] == 0:
# 		rerunids.append(datatable['id'][i])

# rerunids.append(2128515216999497600)
# rerunids.append(2128515496174251008)
# rerunids.append(2075875655448389120)
# ReRunTable['id'] = np.unique(rerunids)

# ReRunTable.write('rerunids.csv', format = 'csv', overwrite = True)
# print(len(ReRunTable))