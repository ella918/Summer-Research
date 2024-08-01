import thejoker as tj
samples = tj.JokerSamples.read('prior_samples_50M.hdf5')
print(samples)

prior = tj.JokerPrior.default( #initializing the default prior
        P_min = 2 * u.day,
        P_max = 1e3 * u.day,
        sigma_K0 = 30 * u.km / u.s,
        sigma_v = 100 * u.km / u.s,
   )
