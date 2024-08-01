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

workpath = '/data2/labs/douglste-laf-lab/mathewea/TheJoker_Outputs' #might have to change this 

unimodal_table = Table.read('/data2/labs/douglste-laf-lab/mathewea/Summer-Research/unimodalcheck.csv')


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

for i in range(len(unimodal_table)):
	if unimodal_table['unimodal'][i] == 1:
		idnum = unimodal_table['id'][i]
		ids.append(idnum)
		if unimodal_table['MCMC'][i] == 0:
			joker_samples = tj.JokerSamples.read(f'{workpath}/{idnum}/rejection_samples_{idnum}.hdf5')
			MCMC.append(0)
		if unimodal_table['MCMC'][i] == 1:
			joker_samples = tj.JokerSamples.read(f'{workpath}/{idnum}/rejection_samples_MCMC_{idnum}.hdf5')
			MCMC.append(1)

		p_median = np.percentile(joker_samples['P'], 50)
		print(p_median)
		p_16 = np.percentile(joker_samples['P'], 16)
		print(p_16)
		p_84 = np.percentile(joker_samples['P'], 84)
		print(p_84)

		e_median1 = np.percentile(joker_samples['e'], 50)
		print(e_median1)
		e_16 = np.percentile(joker_samples['e'], 16)
		print(e_16)
		e_84 = np.percentile(joker_samples['e'], 84)
		print(e_84)

		P_median.append(p_median)
		P_lower.append(p_median - p_16)
		P_upper.append(p_84 - p_median)

		e_median.append(e_median1)
		e_lower.append(e_median - e_16)
		e_upper.append(e_84 - e_median)

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
#data_for_plots['MAP'] = MAP 

data_for_plots.write('data_for_unimodal_plots.csv', format = 'csv', overwrite = True)
print('made table')


#plotting e vs P plot with uncertainties 
fig, ax = plt.subplots()
x = data_for_plots['P_median']
y = data_for_plots['e_median']
#asymmetric_error_x = [data_for_plots['P_lower'], data_for_plots['P_upper']]
#asymmetric_error_y = [data_for_plots['e_lower'], data_for_plots['e_upper']]
#ax.errorbar(x, y, xerr = asymmetric_error_x.all(), yerr = asymmetric_error_y.all())
plt.scatter(x, y)
plt.show()
plt.savefig('Unimodal_Plots')
print('made plot')




