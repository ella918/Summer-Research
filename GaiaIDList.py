from astropy.table import QTable, vstack, unique, Table
import numpy as np

new_6866 = QTable.read('/users/EllaMathews/Research/rcat_ngc6866_v0.fits')
new_6811 = QTable.read('/users/EllaMathews/Research/rcat_ngc6811_v0.fits')
new_ids_6811 = new_6811['GAIAEDR3_ID']
new_ids_6866 = new_6866['GAIAEDR3_ID']

id_list = vstack([new_ids_6811, new_ids_6866])

values, counts = np.unique(id_list, return_counts=True)
id_list_rvs = values[counts >= 3]

idtable = Table()
ids = []
for num in id_list_rvs:
	idnum = np.int64(num)
	ids.append(idnum)
idtable["GAIAEDR3_ID"] = ids

idtable.write("GAIADR3_IDs.csv", format = "csv", overwrite = True)