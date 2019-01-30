# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 11:40:54 2018

@author: nooteboom
"""


import numpy as np
import time
from netCDF4 import Dataset
import math

sp = 's2'#50#'s2'#3
dd = 10
ddeg = 2

dirRead = '/Users/nooteboom/Documents/PhD/parcels/OFES_global/Transtition_matrices/input/'
dirWrite = '/Users/nooteboom/Documents/PhD/parcels/OFES_global/Transtition_matrices/output/box-avgdist/'

#data = Dataset(dirRead + 'timeseries_per_location_inclatlonadv_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd)) + '.nc')
if(type(sp)==str):
    if(sp=='s1'):
        data = Dataset(dirRead + 'timeseries_per_location_inclatlonadv_ddeg%d_spincreaseS1_dd%d'%(ddeg, int(dd)) + '.nc')    
    else:
        data = Dataset(dirRead + 'timeseries_per_location_inclatlonadv_ddeg%d_spincreaseS2_dd%d'%(ddeg, int(dd)) + '.nc')    
else:
    data = Dataset(dirRead + 'timeseries_per_location_inclatlonadv_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd)) + '.nc')
# %%Calculating the transition matrix: 
print 'Calculate the average distance of every grid box'  

lons = data['lon'][:]
lats = data['lat'][:]
lons0 = data['lon0'][:]
lats0 = data['lat0'][:]
vLons = data['vLons'][:]
vLons[vLons==360] -= 360 
vLats = data['vLats'][:]

#%% distance function 
def distance(origin, destination):
    #Calculate distance between two locations in lon en lat
    lat1, lon1 = origin
    lat2, lon2 = destination
    radius = 6371.1 # km

    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c

    return d

#%%
print 'amount of boxes: ', len(vLons)

start = time.time()

nanvalues = 0

trm = np.full((len(vLons)),np.nan) 
for r1 in range(len(vLons)):
    if(r1%1000==0):
        print r1/np.float(len(vLons))
        
    firstts_lon = lons[r1]#[~lons[r1].mask]
    secondts_lon = lons0[r1]#[~lons[r1].mask]
    
    if(len(firstts_lon)>0 and len(secondts_lon)>0):    
        firstts_lat = lats[r1]#[~lats[r1].mask]
        secondts_lat = lats0[r1]#[~lats[r1].mask]
        lon0 = vLons[r1]
        lat0 = vLats[r1]

        dist = np.full((len(firstts_lon)), np.nan)
        
        for j in range(len(firstts_lon)):
            if(firstts_lon[j]>-1000 and firstts_lat[j]>-1000):
#                print 'lon,lats:   ',(lat0,lon0),(firstts_lat[j],firstts_lon[j])
                dist[j] = distance((secondts_lat[j],secondts_lon[j]),(firstts_lat[j],firstts_lon[j]))
#                print dist[j]
        
        trm[r1] = np.nanmean(dist)
 #       if(vLats[r1]<-51 and vLats[r1]>-55):
 #           print vLats[r1]
 #           print trm[r1]
    

    else:
        nanvalues += 1

import matplotlib.pylab as plt

plt.contourf(trm.reshape(len(data['Lats'][:]),len(data['Lons'][:])))
plt.show()

print 'prop. nan values: ', nanvalues / np.float(len(trm))
print 'time \'rTM calculation\' (minutes): ', ((time.time()-start)/60.)
print 'all nan in trm: ', np.isnan(trm).all()

if(type(sp)==str):
    np.savez(dirWrite + 'TM_box-avgdist_ddeg%d_sp'%(ddeg)+sp+'_dd%d'%(int(dd)), TM=np.array(trm), Lons=data['Lons'][:], Lats=data['Lats'][:]) 
else:
    np.savez(dirWrite + 'TM_box-avgdist_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd)), TM=np.array(trm), Lons=data['Lons'][:], Lats=data['Lats'][:]) 