from astropy.table import QTable, vstack

new_6866 = QTable.read('/users/EllaMathews/Research/rcat_ngc6866_v0.fits')
new_6811 = QTable.read('/users/EllaMathews/Research/rcat_ngc6811_v0.fits')
new_ids_6811 = new_6811['GAIAEDR3_ID']
new_ids_6866 = new_6866['GAIAEDR3_ID']

id_list = vstack([new_ids_6811, new_ids_6866])

id_list.write("GAIADR3_IDs.csv", format = 'csv', overwrite = True)