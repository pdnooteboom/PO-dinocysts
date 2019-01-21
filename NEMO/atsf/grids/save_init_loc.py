# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 15:31:22 2017

@author: nooteboom
"""

import numpy as np
from parcels import Field

import time

start = time.time()
for posidx in range(126):
    startid = time.time()

    lonsmin = np.arange(0,359,20) - 180
    lonsmax = np.arange(20,361,20);
    lonsmax = lonsmax - 180
    latsmin = np.array([-75,-50,-25,0, 25, 50, 75])
    latsmax = np.array([-50,-25,0,25, 50, 75, 88])
#Chaning the lon and lat you must also do this within the kernels
    Bind = posidx/18
    minlon = lonsmin[posidx%18]
    maxlon = lonsmax[posidx%18]
    minlat = latsmin[posidx/18]
    maxlat = latsmax[posidx/18]

    sp = 6. #The sinkspeed m/day
    dd = 10. #The dwelling depth

    res = 1 #resolution in degrees

    grid = np.mgrid[int(minlon):int(maxlon):res,int(minlat):int(maxlat):res]+0.5
    n=grid[0].size
    lons = np.reshape(grid[0],n)
    lats = np.reshape(grid[1],n)

#Delete the particles on the land
    filename = 'NEMO-MEDUSA/ORCA0083-N006/domain/bathymetry_ORCA12_V3.3.nc'# Bathymetry file of NEMO
    Land = Field.from_netcdf(filename, variable='Bathymetry',dimensions={'lon':'nav_lon','lat':'nav_lat'})
    Land.grid.check_zonal_periodic()

    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if Land[0,lo,la,0]>dd  ]),np.array([la for lo, la in zip(lons,lats) if Land[0,lo,la,0]>dd ])]

    print 'number of particles beginning: ', len(lats)
#Remove land particles and mediterranean particles
    print 'Time for one posidx (in minutes):', ((time.time()-startid)/60.)

    np.save('coor/releaselocations/lons_id%d_dd%d.npy'%(posidx,int(dd)),lons)
    np.save('coor/releaselocations/lats_id%d_dd%d.npy'%(posidx,int(dd)),lats)

print 'Time of the whole execution (in minutes):', ((time.time()-start)/60.)
