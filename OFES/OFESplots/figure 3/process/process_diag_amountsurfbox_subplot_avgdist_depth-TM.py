# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 15:43:10 2019

@author: nooteboom
"""
import numpy as np
from parcels import Field
import math
from netCDF4 import Dataset

def GetOFESLandArray(filename, fieldname):
    """
    Function to return a Field with 1's at land points and 0's at ocean points, based on python basemap
    :param f: a field of the .nc file, but not a parcels field object! For OFES, land points are masked. This is used here!
    """
    pfile = Dataset(filename, 'r')
    Lon = pfile.variables['LONN1799_1800'][:]
    Lat = pfile.variables['LAT'][:]
    f = pfile.variables[fieldname][:]
    f = f[0,0,:,:] 
    Land=Field('Land',f,transpose=False, lon=Lon,lat=Lat)
    return Land

#Interpolates the bathymetry:
Bath = GetOFESLandArray('topography_OFES.nc', 'HT')

ddeg = 2 # resolution of the binning
sp = 6#'s2'
dd = 100
res = 1

tmdir = '/OFES/OFESres/TM/'

if(type(sp)==str):
    data = np.load(tmdir + 'box-box/TMglobal_bin'+str(ddeg)+'_dd'+str(int(dd)) +'_sp'+sp+"_res"+str(res) + '.npz')
else:
    data = np.load(tmdir + 'box-box/TMglobal_bin'+str(ddeg)+'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res) + '.npz')
TM = data['TM'][:]
Lons = data['Lons'][:]
Lats = data['Lats'][:]
#%%Load avg drift distance
if(type(sp)==str):
    datadist = np.load(tmdir + 'box-avgdist/TM_box-avgdist_ddeg%d_sp'%(ddeg)+sp+'_dd%d.npz'%(dd))
else:
    datadist = np.load(tmdir + 'box-avgdist/TM_box-avgdist_ddeg%d_sp%d_dd%d.npz'%(ddeg, sp, dd))

TMd = datadist['TM'][:]
    
Lonsd = datadist['Lons'][:]
Latsd = datadist['Lats'][:]
vLons, vLats = np.meshgrid(Lonsd, Latsd)
vLons = vLons.flatten(); vLats = vLats.flatten();

Lonsd[Lonsd>=360] -= 360; 

TMd = TMd[np.logical_and(vLons<360,vLats<=73)]
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
        lo = x[i]+1
        if(lo>=180):
            lo -= 360
            if(lo>180):
                lo -= 180
        if(~np.ma.is_masked(Bath[0,lo,y[i]+1,0])):
            depths[i] = Bath[0,lo,y[i]+1,0]    
                
    return depths

vLons, vLats = np.meshgrid(Lons, Lats)
vLons = vLons.flatten(); vLats = vLats.flatten();

surface = np.full(vLons.shape[0], 0)

for i in range(len(vLons)):
    lonmin = vLons[i]; lonmax = vLons[i]+ddeg;
    latmin = vLats[i]; latmax = vLats[i]+ddeg;   
    dis1 = distance((latmin,lonmin),(latmax,lonmin))
    dis2 = distance((latmin,lonmin),(latmin, lonmax))
    surface[i] = dis1 * dis2
        
Lons = data['Lons'][:]
Lats = data['Lats'][:]
vLons, vLats = np.meshgrid(Lons, Lats)
vLons = vLons.flatten(); vLats = vLats.flatten();
depth = box_depths(vLons, vLats)
        
threshold = 0.
connected = (TM>threshold).astype(int)

for i in range(connected.shape[0]):
    connected[i] = connected[i] * surface

            
np.savez('toplot_TM_diag_amountboxes_avgdist_scatter_sp'+str(sp)+'_dd'+str(dd), 
         TM = np.diagonal(TM, offset=0), connected = np.sum(connected, axis=0), 
        TMd = TMd, Lons=Lons, Lats=Lats, 
        depth = depth)           
            