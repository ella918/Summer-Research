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

DATA_PATH = os.getenv("DATA_PATH", "/users/EllaMathews/Summer-Research") #environment variable 
prior = tj.JokerSamples.read('/data2/labs/douglste-laf-lab/mathewea/TheJoker_Outputs/2128128296283380480/prior_samples_10.0M_2128128296283380480.hdf5')
rejectionsamples = tj.JokerSamples.read('/data2/labs/douglste-laf-lab/mathewea/TheJoker_Outputs/2128128296283380480/rejection_samples_2128128296283380480.hdf5')
new_6866 = QTable.read('/data2/labs/douglste-laf-lab/mathewea/rcat_ngc6866_v0.fits')
new_6811 = QTable.read('/data2/labs/douglste-laf-lab/mathewea/rcat_ngc6811_v0.fits')
id_list = Table.read('/data2/labs/douglste-laf-lab/mathewea/Summer-Research/GAIADR3_IDs.csv')

id_num = 2128128296283380480
new_ids_6811 = new_6811['GAIAEDR3_ID']
new_ids_6866 = new_6866['GAIAEDR3_ID']
datamatched6811 = new_6811[id_num == new_ids_6811]
datamatched6866 = new_6866[id_num == new_ids_6866]
matched = vstack([datamatched6811, datamatched6866])
t1 = Time(matched["DATE-OBS"], format = "fits", scale = "tcb")
data = tj.RVData(t = t1, rv = matched['vrad']*(u.kilometer/u.second), rv_err = matched['vrad_err']*(u.kilometer/u.second)) 

plt.rcParams.update({'font.size': 14})

fig1, ax1 =plt.subplots()
_ = data.plot()
ax1.set_title("Inital RV Data")
fig1.savefig('initialRV_AAS.png')
print('RV data plotted')

fig2, ax2 = plt.subplots()
ax2.scatter(prior["P"], prior["e"], s=20, lw=0, alpha=0.5)
ax2.set_xscale("log")
ax2.set_xlim(1, 1e3)
ax2.set_ylim(0, 1)
ax2.set_xlabel("$P$ [day]")
ax2.set_ylabel("$e$")
ax2.set_title("10 Million Prior Samples")
fig2.savefig('priors_AAS.png')
print('Prior Plotted')


fig3, ax3 = plt.subplots()
with quantity_support():
    ax3.scatter(rejectionsamples["P"], rejectionsamples["e"], s=20, lw=0, alpha=0.5)
    ax3.set_xscale("log")
    ax3.set_xlim(1, 1e3)
    ax3.set_ylim(0, 1)
    ax3.set_xlabel("$P$ [day]")
    ax3.set_ylabel("$e$")
    ax3.set_title("Rejection Samples")
fig3.savefig('rejectionsamples_AAS.png')
print('rejection samples Plotted')

fig4, ax4 = plt.subplots()
_ = tj.plot_rv_curves(rejectionsamples, data=data) 
ax4.set_title("Possible Orbits")
fig4.savefig('posorbits_AAS.png')
print("RV curves after rejection plotted")





