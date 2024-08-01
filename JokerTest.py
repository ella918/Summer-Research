import thejoker as tj
import astropy.table as at
import astropy.units as u
import matplotlib.pyplot as plt
# import numpy as np
# import h5py
# import thejoker as tj
# from astropy.table import QTable, Table, Column, vstack, unique
# from astropy.time import Time
# from astropy.visualization.units import quantity_support
# import astropy.coordinates as coord
# import pymc as pm
# import arviz as az
# import argparse
# import os
# import schwimmbad


samples = tj.JokerSamples.read('prior_samples_50M.hdf5')
print(samples)

prior = tj.JokerPrior.default( #initializing the default prior
        P_min = 2 * u.day,
        P_max = 1e3 * u.day,
        sigma_K0 = 30 * u.km / u.s,
        sigma_v = 100 * u.km / u.s,
   )
