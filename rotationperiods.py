#imports 
import astropy.io.ascii as at
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import h5py
from astropy.table import QTable, Table, Column, join, unique, vstack
from astropy.time import Time
from astropy.visualization.units import quantity_support
from astropy import table
import os

#importing catalogs
rot6811 = at.read("/users/EllaMathews/Summer-Research/ngc6811_rotation_curtis2019.txt")
rot6866 = at.read("/users/EllaMathews/Summer-Research/ngc6866_rotation_balona2013.tsv",data_start=3,delimiter="|")
rot_long = at.read("/users/EllaMathews/Summer-Research/clusters_rotation_long2023.tsv",data_start=3,delimiter="|")

data_6811 = QTable.read('/users/EllaMathews/Research/rcat_ngc6811_v0.fits')
data_6866 = QTable.read('/users/EllaMathews/Research/rcat_ngc6866_v0.fits')
old_6866 = QTable.read('/users/EllaMathews/Summer-Research/Combined_Data_6866.csv')
old_6811 = QTable.read('/users/EllaMathews/Research/Combined_Data_6811.csv')
values6811, counts6811 = np.unique(data_6811['GAIAEDR3_ID'], return_counts = True)
values6866, counts6866 = np.unique(data_6866['GAIAEDR3_ID'], return_counts = True)

rotlong_newnames = Table()
rotlong_newnames['Gaia'] = rot_long['Gaia']
rotlong_newnames['Per'] = rot_long['PRot']
rotlong_newnames['e_Per'] = rot_long['e_PRot']

all6811 = Table()
all6811['Gaia'] = values6811
all6811['num_RV'] = counts6811

# combinedrot6811 = vstack([rot6811, rotlong_newnames])
# combinedrot6811_P = combinedrot6811['Per']
# combinedrot6811_withP = combinedrot6811[combinedrot6811_P > 0]
# rotRV6811 = join(combinedrot6811, all6811, keys = 'Gaia')
# rotRV6811_wantedcolumns = rotRV6811['Gaia','Per', 'e_Per', 'num_RV']
# print(rotRV6811_wantedcolumns)
# print(data_6811.info)

old_gaia_6811 = []
old_KIC_6811 = []
for name in old_6811['DR3Name']:
	num = np.int64(name[9:])
	old_gaia_6811.append(num)
for idnum in old_6811['id']:
	old_KIC_6811.append(idnum)
KICGaia6811 = Table()
KICGaia6811['KIC'] = old_KIC_6811
KICGaia6811['Gaia'] = old_gaia_6811
bothids6811 = table.unique(KICGaia6811, keys = 'Gaia')

rot6811_wantedcolumns = rot6811['KIC', 'Per', 'e_Per']
rot6811bothids = join(rot6811_wantedcolumns, bothids6811, keys ='KIC')
rotcombined6811 = vstack([rot6811bothids, rotlong_newnames])
rotcombined6811_P = rotcombined6811['Per']
rotcombined6811_withP = rotcombined6811[rotcombined6811_P > 0]

rotRV6811 = join(rotcombined6811_withP, all6811, keys = 'Gaia')

colordata6811 = Table()
colordata6811['Gaia'] = data_6811['GAIAEDR3_ID']
colordata6811['BP-RP'] = data_6811['GAIAEDR3_BP'] - data_6811['GAIAEDR3_RP']
colordata6811['Gmag'] = data_6811['GAIAEDR3_G']
colordata6811['e_BP-RP'] = data_6811['GAIAEDR3_BP_ERR'] - data_6811['GAIAEDR3_RP_ERR']
colordata6811['e_Gmag'] = data_6811['GAIAEDR3_G_ERR']

combinedwithcolor6811 = table.unique(join(rotRV6811, colordata6811, keys = 'Gaia'), keys = 'Gaia')
wantedcolumns6811 = combinedwithcolor6811['Gaia', 'Per', 'BP-RP', 'Gmag', 'e_BP-RP', 'e_Gmag', 'num_RV']
wantedcolumns6811.write('RotationPeriodRV_6811.csv', format = 'csv', overwrite = True)

all6866 = Table()
all6866['Gaia'] = values6866
all6866['num_RV'] = counts6866

rotRV6866 = join(rotlong_newnames, all6866, keys = 'Gaia')

colordata6866 = Table()
colordata6866['Gaia'] = data_6866['GAIAEDR3_ID']
colordata6866['BP-RP'] = data_6866['GAIAEDR3_BP'] - data_6866['GAIAEDR3_RP']
colordata6866['Gmag'] = data_6866['GAIAEDR3_G']
colordata6866['e_BP-RP'] = data_6866['GAIAEDR3_BP_ERR'] - data_6866['GAIAEDR3_RP_ERR']
colordata6866['e_Gmag'] = data_6866['GAIAEDR3_G_ERR']

combinedwithcolor6866 = table.unique(join(rotRV6866, colordata6866, keys = 'Gaia'), keys = 'Gaia')
wantedcolumns6866 = combinedwithcolor6866['Gaia', 'Per', 'BP-RP', 'Gmag', 'e_BP-RP', 'e_Gmag', 'num_RV']
wantedcolumns6866.write('RotationPeriodRV_6866.csv', format = 'csv', overwrite = True)

counterPer6811 = 0 	
counterPer3RV6811 = 0
counterPer8RV6811 = 0

counterPer6866 = 0
counterPer3RV6866 = 0
counterPer8RV6866 = 0

counterG6811 = 0
counterG3RV6811 = 0
counterG8RV6811 = 0

counterG6866 = 0
counterG3RV6866 = 0
counterG8RV6866 = 0

for i in range(len(wantedcolumns6811)):
	if wantedcolumns6811['BP-RP'][i] >= 0.4 and wantedcolumns6811['BP-RP'][i] <= 1.4:
		if wantedcolumns6811['Per'][i] > 0 and wantedcolumns6811['Per'][i] <=20:
			counterPer6811 = counterPer6811+1
			if wantedcolumns6811['num_RV'][i] >= 3:
				counterPer3RV6811 = counterPer3RV6811 + 1
			if wantedcolumns6811['num_RV'][i] >= 8:
				counterPer8RV6811 = counterPer8RV6811 + 1
print('counterPer6811 = ', counterPer6811)
print('counterPer3RV6811 = ', counterPer3RV6811)
print('counterPer8RV6811 = ', counterPer8RV6811)
print('\n')

for i in range(len(wantedcolumns6811)):
	if wantedcolumns6811['BP-RP'][i] >= 0.4 and wantedcolumns6811['BP-RP'][i] <= 1.4:
		if wantedcolumns6811['Gmag'][i] > 12 and wantedcolumns6811['Gmag'][i] <=17:
			counterG6811 = counterG6811+1
			if wantedcolumns6811['num_RV'][i] >= 3:
				counterG3RV6811 = counterG3RV6811 + 1
			if wantedcolumns6811['num_RV'][i] >= 8:
				counterG8RV6811 = counterG8RV6811 + 1
print('counterG6811 = ', counterG6811)
print('counterG3RV6811 = ', counterG3RV6811)
print('counterG8RV6811 = ', counterG8RV6811)
print('\n')

for i in range(len(wantedcolumns6866)):
	if wantedcolumns6866['BP-RP'][i] >= 0.4 and wantedcolumns6866['BP-RP'][i] <= 1.4:
		if wantedcolumns6866['Per'][i] > 0 and wantedcolumns6866['Per'][i] <=20:
			counterPer6866 = counterPer6866+1
			if wantedcolumns6866['num_RV'][i] >= 3:
				counterPer3RV6866 = counterPer3RV6866 + 1
			if wantedcolumns6866['num_RV'][i] >= 8:
				counterPer8RV6866 = counterPer8RV6866 + 1
print('counterPer6866 = ', counterPer6866)
print('counterPer3RV6866 = ', counterPer3RV6866)
print('counterPer8RV6866 = ', counterPer8RV6866)
print('\n')

for i in range(len(wantedcolumns6866)):
	if wantedcolumns6866['BP-RP'][i] >= 0.4 and wantedcolumns6866['BP-RP'][i] <= 1.4:
		if wantedcolumns6866['Gmag'][i] > 12 and wantedcolumns6866['Gmag'][i] <=17:
			counterG6866 = counterG6866+1
			if wantedcolumns6866['num_RV'][i] >= 3:
				counterG3RV6866 = counterG3RV6866 + 1
			if wantedcolumns6866['num_RV'][i] >= 8:
				counterG8RV6866 = counterG8RV6866 + 1
print('counterG6866 = ', counterG6866)
print('counterG3RV6866 = ', counterG3RV6866)
print('counterG8RV6866 = ', counterG8RV6866)
print('\n')




















fig, axes = plt.subplots(nrows = 2, ncols = 2, sharex = True, sharey = 'row', width_ratios = [1, 1.25])

ax1 = axes[0,0]
x1 = wantedcolumns6811['BP-RP']
y1 = wantedcolumns6811['Per']
scatter1=ax1.scatter(x1, y1, c = wantedcolumns6811['num_RV'], marker = 'o', alpha = 0.5, vmin = 1, vmax = 10, s =10)
ax1.set_xlim(0.4, 1.4)
ax1.set_ylim(0, 20)
ax1.set_title('NGC 6811')
ax1.tick_params(axis = 'x', labelbottom=False)
ax1.set_ylabel('Rotation Period (d)')
#plt.colorbar(scatter1)

ax2 = axes[1,0]
x2 = wantedcolumns6811['BP-RP']
y2 = wantedcolumns6811['Gmag']
scatter2=ax2.scatter(x2, y2, c = wantedcolumns6811['num_RV'], marker = 'o', alpha = 0.5, vmin =1, vmax=10, s =10)
ax2.set_ylim(12, 17)
ax2.set_xlabel('BP-RP')
ax2.set_ylabel('G')
ax2.invert_yaxis()
#plt.colorbar(scatter2)


ax3 = axes[0,1]
x3 = wantedcolumns6866['BP-RP']
y3 = wantedcolumns6866['Per']
scatter3=ax3.scatter(x3, y3, c = wantedcolumns6866['num_RV'], marker = 'o', alpha = 0.5, vmin=1, vmax=10, s =10)
ax3.set_title('NGC 6866')
plt.colorbar(scatter3)

ax4 = axes[1,1]
x4 = wantedcolumns6866['BP-RP']
y4 = wantedcolumns6866['Gmag']
scatter4=ax4.scatter(x4, y4, c = wantedcolumns6866['num_RV'], marker = 'o', alpha = 0.5, vmin=1, vmax=10, s =10)
ax4.set_xlabel('BP-RP')

plt.colorbar(scatter4)

plt.subplots_adjust(hspace = .1, wspace = .105)
#plt.show()
#plt.savefig('NewPlots', bbox_inches='tight')

		
