# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 11:04:41 2018

Script which calculates the location of the ACC.
Within the OFES data, the maximum flow speed for every latitude is considered 
within some latitude range. Then a running mean is computed to smooth the ACC
trajectory.

@author: nooteboom
"""

from netCDF4 import Dataset
import numpy as np
from os import path
from datetime import date
from datetime import timedelta as delta


def snapshot_function(start, end, delta):
    """
    
    """
    curr = start
    result = []
    result.append('{:%Y%m%d}'.format(curr))
    while curr <= end:
     #   yield curr
        curr += delta
        result.append('{:%Y%m%d}'.format(curr))
    del result[-1]
    result = np.array(result)
    return result

directory = '/projects/0/palaeo-parcels/OFESdata/'
snapshots = snapshot_function(date(2000, 1, 3), date(2005, 12, 29), delta(days=3))

ufiles = [path.join(directory, 'y'+s[:4], 'u_vel', "nest_1_"+s+"000000u.nc".format(s)) for s in snapshots]
vfiles = [path.join(directory, 'y'+s[:4], 'v_vel', "nest_1_"+s+"000000v.nc".format(s)) for s in snapshots]

tfiles = [path.join(directory, 'y'+s[:4], 'temp', "nest_1_"+s+"000000t.nc".format(s)) for s in snapshots]

maxlat = 380
minlat = 100

for i in range(len(ufiles)):
    if(i%10==0):
        print i/float(len(ufiles))
    uf = Dataset(ufiles[i])
    vf = Dataset(vfiles[i])
    
    tf = Dataset(tfiles[i])

    if(i==0):
        u = uf['zu'][0,0,minlat:maxlat,:]
        u[np.isnan(u)] = 0
        v = vf['zv'][0,0,minlat:maxlat,:]
        v[np.isnan(v)] = 0
        print u.shape
        print uf['Latitude'][minlat:maxlat]
        
        t = tf['temperature'][0,0,minlat:maxlat,:]
    else:
        ut = uf['zu'][0,0,minlat:maxlat,:]
        ut[np.isnan(ut)] = 0
        ut[ut>10000] = 0
        u += ut#uf['zu'][0,0,minlat:maxlat,:]
        vt = vf['zv'][0,0,minlat:maxlat,:]
        vt[np.isnan(vt)] = 0
        vt[vt>10000] = 0
        v += vt

        tt = tf['temperature'][0,0,minlat:maxlat,:]
        tt[np.isnan(tt)] = 0
       
        t+= tt
        
        
#Calculate average velocity fields: 
u = u / np.float(len(ufiles))
v = v / np.float(len(vfiles))

t = t / np.float(len(tfiles))

#Calculate magnitude of the velocity
ub = np.sqrt(np.power(u,2) + np.power(v,2))

acclat = np.full(ub.shape[1], np.nan)
acclon = np.full(ub.shape[1], np.nan)

SAFlat = np.full(ub.shape[1], np.nan)
SAFlon = np.full(ub.shape[1], np.nan)

for loni in range(ub.shape[1]):
    acclon[loni] = uf['Longitude'][loni]
    lon = uf['Longitude'][loni]    
    lati = np.argmax(u[:,loni])
    acclat[loni] = uf['Latitude'][minlat+lati]
    
    SAFlon[loni] = tf['Longitude'][loni]
    lon = tf['Longitude'][loni]    
    lati = np.argmin((np.abs(t[:,loni])-8))
    SAFlat[loni] = tf['Latitude'][minlat+lati]    
    
np.savez('ACClocation.npz', acclon = acclon, acclat = acclat)
np.savez('SAFlocation.npz', SAFlon = SAFlon, SAFlat = SAFlat)
    