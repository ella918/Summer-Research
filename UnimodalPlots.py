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

workpath = '/data2/labs/douglste-laf-lab/mathewea/200.0M_new'
unimodal_table = Table.read('/data2/labs/douglste-laf-lab/mathewea/Summer-Research/unimodalcheck_200M_new.csv')

data_for_plots = QTable()
ids = []
P_median = []
P_lower = []
P_upper = []
P_std = []
e_median = []
e_lower = []
e_upper = []
e_std = []
MCMC = []
uni =  []

for i in range(len(unimodal_table)):
	if unimodal_table['unimodal'][i] == 1 or unimodal_table['unimodal'][i] == 0:
		idnum = unimodal_table['id'][i]
		ids.append(idnum)
		if unimodal_table['unimodal'][i] == 1:
			uni.append(1)
		if unimodal_table['unimodal'][i] == 0:
			uni.append(0)
		if unimodal_table['MCMC'][i] == 0:
			joker_samples = tj.JokerSamples.read(f'{workpath}/{idnum}/rejection_samples_200.0M_{idnum}_new.hdf5')
			MCMC.append(0)
		if unimodal_table['MCMC'][i] == 1:
			joker_samples = tj.JokerSamples.read(f'{workpath}/{idnum}/rejection_samples_MCMC_200.0M_{idnum}_new.hdf5')
			MCMC.append(1)

		p_median = np.percentile(joker_samples['P'], 50)
		#print(p_median)
		p_16 = np.percentile(joker_samples['P'], 16)
		#print(p_16)
		p_84 = np.percentile(joker_samples['P'], 84)
		#print(p_84)
		p_std = np.std(joker_samples['P'])

		e_median1 = np.percentile(joker_samples['e'], 50)
		#print(e_median1)
		e_16 = np.percentile(joker_samples['e'], 16)
		#print(e_16)
		e_84 = np.percentile(joker_samples['e'], 84)
		#print(e_84)
		e_std1 = np.std(joker_samples['e'])

		P_median.append(p_median)
		P_lower.append(p_median - p_16)
		P_upper.append(p_84 - p_median)
		P_std.append(p_std)

		e_median.append(e_median1)
		e_lower.append(e_median1 - e_16)
		e_upper.append(e_84 - e_median1)
		e_std.append(e_std1)
		#Map = tj.MAP_sample(joker_samples)
		#MAP.append(Map)

data_for_plots['id'] = ids
data_for_plots['MCMC'] = MCMC
data_for_plots['unimodal']= uni
data_for_plots['P_median'] = P_median
data_for_plots['P_lower'] = P_lower
data_for_plots['P_upper'] = P_upper 
data_for_plots['P_std'] = P_std
data_for_plots['e_median'] = e_median
data_for_plots['e_lower'] = e_lower
data_for_plots['e_upper'] = e_upper 
data_for_plots['e_std'] = e_std
data_for_plots.write('data_for_plots_all_200M_new.csv', format = 'csv', overwrite = True)
print('made table')

data_for_plots['med/std'] = data_for_plots['P_median']/data_for_plots['P_std']
check = data_for_plots['med/std'] >3
x = data_for_plots['P_median'].value
y = data_for_plots['e_median'].value
xerr_low = data_for_plots['P_lower'].value
xerr_up = data_for_plots['P_upper'].value
yerr_low = data_for_plots['e_lower'].value
yerr_up = data_for_plots['e_upper'].value

print(type(x),type(y),type(xerr_low))
print(type(xerr_up),type(yerr_low),type(yerr_up))

#plotting e vs P plot with uncertainties 
fig, ax = plt.subplots()
asymmetric_error_x = np.array([np.array(xerr_low[check]), np.array(xerr_up[check])])
#print(asymmetric_error_x)
asymmetric_error_y = np.array([np.array(yerr_low[check]), np.array(yerr_up[check])])
#print(asymmetric_error_y)
ax.errorbar(x[check], y[check], xerr = asymmetric_error_x, yerr = asymmetric_error_y, fmt = 'o')
ax.set_xlabel('P (d)')
ax.set_xscale('log')
ax.set_ylabel('e')
ax.set_title('e vs P for Stars P/std >3  in ngc6811 and ngc6866')
plt.show()
plt.savefig('evP_check_std_200M_new')
print('made plot')

for i in range(len(data_for_plots)):
   if np.float64(data_for_plots['P_median'][i]) < 20:
      print(data_for_plots['id'][i])
      print(data_for_plots['P_median'][i])
      print(data_for_plots['e_median'][i])


