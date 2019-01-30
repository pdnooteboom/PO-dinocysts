# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 11:40:54 2018

@author: nooteboom
"""


import numpy as np
import time
from netCDF4 import Dataset

sp = 500#'s2'#250
dd = 10
ddeg = 1

n = 5

dirRead = '/Users/nooteboom/Documents/PhD/parcels/OFES_global/Transtition_matrices/input/'#'/projects/0/palaeo-parcels/OFESres/particledata/sp%d_dd%d/'%(int(sp),int(dd))
dirWrite = '/Users/nooteboom/Documents/PhD/parcels/OFES_global/Transtition_matrices/output/box-percentile/'

var = 'temp'   #Put to 'temp' or 'salin'

if(type(sp)==str):
    if(sp=='s1'):
        data = Dataset(dirRead + 'timeseries_per_location_inclatlonadv_ddeg%d_spincreaseS1_dd%d'%(ddeg, int(dd)) + '.nc')    
    else:
        data = Dataset(dirRead + 'timeseries_per_location_inclatlonadv_ddeg%d_spincreaseS2_dd%d'%(ddeg, int(dd)) + '.nc')    
else:
    data = Dataset(dirRead + 'timeseries_per_location_inclatlonadv_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd)) + '.nc')

# %%Calculating the transition matrix: 
print 'Calculate the wasserstein distance of every grid box'  

ts = data[var][:]
tsfix = data['fix'+var][:]
vLons = data['vLons'][:]

#%%
print 'amount of boxes: ', len(vLons)

start = time.time()

nanvalues = 0
print ts.shape, tsfix.shape
trm = np.full((len(vLons)),np.nan) 
for r1 in range(len(vLons)):
    if(r1%1000==0):
        print r1/np.float(len(vLons))
    firstts = ts[r1][~ts[r1].mask]
    secondts = tsfix[r1][~tsfix[r1].mask]
    if(len(firstts)>0 and len(secondts)>0):
        trm[r1] = np.nanpercentile(firstts,n) - np.nanpercentile(secondts,n)
    else:
        nanvalues += 1

print 'prop. nan values: ', nanvalues / np.float(len(trm))
print 'time \'rTM calculation\' (minutes): ', ((time.time()-start)/60.)
print 'all nan in trm: ', np.isnan(trm).all()

if(type(sp)==str):
    np.savez(dirWrite + 'TM_box-%dpercentiledif_'%(n)+var+'_ddeg%d_sp'%(ddeg)+sp+'_dd%d'%(int(dd)), TM=np.array(trm), Lons=data['Lons'][:], Lats=data['Lats'][:]) 
else:
    np.savez(dirWrite + 'TM_box-%dpercentiledif_'%(n)+var+'_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd)), TM=np.array(trm), Lons=data['Lons'][:], Lats=data['Lats'][:]) 
