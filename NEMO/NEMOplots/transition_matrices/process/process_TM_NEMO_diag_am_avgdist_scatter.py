# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 14:49:22 2019

@author: nooteboom
"""





import numpy as np
import matplotlib.pylab as plt
import matplotlib
from netCDF4 import Dataset
import math

from parcels import Field

filename = 'bathymetry_ORCA12_V3.3.nc'
Bath = Field.from_netcdf(filename, variable='Bathymetry',dimensions={'lon':'nav_lon','lat':'nav_lat'})
#Bath.grid.check_zonal_periodic()

# True if surface instead of amount of surface grid boxes
surfacebool = True

ddeg = 2 # resolution of the binning
sp = 6
dd = 10
res = 1

tmdir = '/Users/nooteboom/Documents/PhD/parcels/NEMO/atsf/Transition_matrices/'

data = np.load(tmdir + 'output/box-box/TMglobal_bin'+str(ddeg)+'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res) + '.npz')
TM = data['TM'][:]
Lons = data['Lons'][:]
Lats = data['Lats'][:]
#%%Load avg drift distance
datadist = np.load(tmdir + 'output/box-avgdist/TM_box-avgdist_ddeg%d_sp%d_dd%d.npz'%(ddeg, sp, dd))
TMd = datadist['TM'][:]
Lonsd = datadist['Lons'][:]
Latsd = datadist['Lats'][:]
vLons, vLats = np.meshgrid(Lonsd, Latsd)
vLons = vLons.flatten(); vLats = vLats.flatten();

#Lonsd[Lonsd>=180] -= 360; 

#TMd = TMd[np.logical_and(vLons>-0.5,vLats<=86)]
#%% Calculate array that gives every grid box a surface value
def distance(origin, destination):
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

def box_depths(x, y):
    depths = np.full(x.shape, np.nan)
    
    for i in range(x.shape[0]):
        if i%1000==0:
            print i/np.float(x.shape[0])
#                print i,j
        lo = x[i]+1
#        print 'lo  ',lo
        if(lo>=180):
            mul = int(lo/180)
            mul = int(np.ceil(mul/2.))
#                    print lo
#                    print mul
            lo -= mul*360
        elif(lo<=-180):
            mul = int(-lo/180)
            mul = int(np.ceil(mul/2.))
#                    print lo
#                    print mul                    
            lo += mul*360                    
        if(Bath[0,lo,y[i]+1,0]!=0):
            depths[i] = Bath[0,lo,y[i]+1,0]                    
    return depths

if(surfacebool):
    vLons, vLats = np.meshgrid(Lons, Lats)
    vLons = vLons.flatten(); vLats = vLats.flatten();

    surface = np.full(vLons.shape[0], 0)
    
    for i in range(len(vLons)):
        lonmin = vLons[i]; lonmax = vLons[i]+ddeg;
        latmin = vLats[i]; latmax = vLats[i]+ddeg;   
        dis1 = distance((latmin,lonmin),(latmax,lonmin))
        dis2 = distance((latmin,lonmin),(latmin, lonmax))
        surface[i] = dis1 * dis2


threshold = 0.
connected = (TM>threshold).astype(int)   
for i in range(connected.shape[0]):
    connected[i] = connected[i] * surface

vLons, vLats = np.meshgrid(Lonsd, Latsd)
vLons = vLons.flatten(); vLats = vLats.flatten()
depthd = box_depths(vLons, vLats)

vLons, vLats = np.meshgrid(Lons, Lats)
vLons = vLons.flatten(); vLats = vLats.flatten()
depth = box_depths(vLons, vLats)



print TM.shape, connected.shape, depth.shape
print TM.shape, TMd.shape, depthd.shape

np.savez('toplot_TM_NEMO_diag_amountboxes_avgdist_scatter', 
         TM = np.diagonal(TM, offset=0), connected = np.sum(connected, axis=0), 
        TMd = TMd, Lons=Lons, Lats=Lats, 
        Lons2=Lonsd, Lats2=Latsd, depth = depth, depthd = depthd)