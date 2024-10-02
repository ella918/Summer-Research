
from astropy.table import QTable, vstack 
import os

id_num = 2128124963389008384
DATA_PATH = os.getenv('DATA_PATH')
workpath = f'/{DATA_PATH}/'
new_6866 = QTable.read(f'{workpath}/rcat_ngc6866_v0.fits')
new_6811 = QTable.read(f'{workpath}/rcat_ngc6811_v0.fits')
new_ids_6811 = new_6811['GAIAEDR3_ID']
new_ids_6866 = new_6866['GAIAEDR3_ID']
datamatched6811 = new_6811[id_num == new_ids_6811]
datamatched6866 = new_6866[id_num == new_ids_6866]
matched = vstack([datamatched6811, datamatched6866])

matched = matched["DATE-OBS","vrad","vrad_err"]

matched.write('MWE_RVData.csv')