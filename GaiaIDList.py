from astropy.table import QTable, vstack, unique, Table
import numpy as np

new_6866 = QTable.read('/users/EllaMathews/Research/rcat_ngc6866_v0.fits')
new_6811 = QTable.read('/users/EllaMathews/Research/rcat_ngc6811_v0.fits')
new_ids_6811 = new_6811['GAIAEDR3_ID']
new_ids_6866 = new_6866['GAIAEDR3_ID']

id_list = vstack([new_ids_6811, new_ids_6866])

values, counts = np.unique(id_list, return_counts=True)
id_list_rvs = values[counts >= 3]

id_list_empty_rvs = [2080063763947694976, 2080090255310324736, 2080091492261057408, 2081881875143510016, 2128099880779629056, 2128108294615662720, 2128112452142469376, 2128126539637279360, 2128131319942169088, 2128131384359936768, 2128132002840746752, 2128141387345637504, 2128158154896484480, 2128170180806442880, 2128170764920170624, 2128198076613025408]

idtable = Table()
ids = []
for num in id_list_rvs:
	idnum = np.int64(num)
	if idnum not in id_list_empty_rvs:
		ids.append(idnum)
idtable["GAIAEDR3_ID"] = ids
print(len(idtable))
idtable.write("GAIADR3_IDs.csv", format = "csv", overwrite = True)