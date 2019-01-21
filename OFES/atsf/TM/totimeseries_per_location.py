# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 16:27:01 2018

@author: nooteboom
"""
import numpy as np
from numba import jit
import time
from netCDF4 import Dataset

pltse = 'annual'

def find_down(array,value):
    if(value>=array[0]):
        array = [n-value for n in array]
        idx = np.array([n for n in array if n<0]).argmax()
    else:
        idx = np.nan
    return idx  

#%% To open a file
res = 1
sp = 500#'s2'#3.
dd = 10. 

print 'sp:  ', sp

ddeg = 1 # binning resolution of the bottom boxes

#dirRead is where to read the data from. dirWrite is where to write the data
if(type(sp)==str):
    if(sp=='s1'):
        dirRead = '/OFES/OFESres/particledata/spincrease_scenario/'
        dirWrite = 'timeseries_per_location/ncfiles/spincrease_scenario/'
    else:
        dirRead = '/OFES/OFESres/particledata/spincrease_scenario2/'
        dirWrite = 'timeseries_per_location/ncfiles/spincrease_scenario2/'
else:
    dirRead = '/OFES/OFESres/particledata/sp%d_dd%d/'%(int(sp),int(dd))
    dirWrite = 'timeseries_per_location/ncfiles/sp%d_dd%d/'%(int(sp),int(dd))

dirRead2 = '/OFES/OFESres/particledata/surface/'

adv = True
nadv = True

#%%
if(type(sp)==str):
    if(sp=='s1'):
        pfile = Dataset(dirRead + 'concatenated_spS1_dd%d_res%d.nc'%(int(dd),res) )
    else:
        pfile = Dataset(dirRead + 'concatenated_sps2_dd%d_res%d.nc'%(int(dd),res) )
else:
    pfile = Dataset(dirRead + 'concatenated_sp%d_dd%d_res%d.nc'%(int(sp),int(dd),res) )

pfilefix = Dataset(dirRead2 + 'concatenatedsurface_dd%d_res%d.nc'%(int(dd),res) )

#%%
minlon = min(pfile['lon0'][:])
maxlon = max(pfile['lon0'][:])+1
minlat = min(pfile['lat0'][:])
maxlat = max(pfile['lat0'][:])+1
#%% Now loop over all regions to define the T matrix
    
#Define the grid at the bottom:
Lons = np.arange(minlon,maxlon+ddeg,ddeg)-0.5
Lats = np.arange(minlat,maxlat+ddeg,ddeg)-0.5
Lonss, Latss = np.meshgrid(Lons, Lats)

vLons = Lonss.reshape(Lonss.size) 
vLats = Latss.reshape(Latss.size)   

@jit
def makets_ex():
    ts = []
    for i in range(len(vLons)): ts.append([]);
    return ts
    
ts = makets_ex()

print len(ts)

def construct_ts(tstemp, tssalin, tslon, tslat, tslon0, tslat0, tsage, lats0, lons0,  vLons, vLats, Lons, Lats, temp, salin, lonadv, latadv, ageadv):

    """Return a list with first dimension len(VLons) for all the bottom grid boxes
    and second dimension the temperatures at the surface related to the bottom box"""  
    le = len(lons0)
    
    maxlents = 0    
    #Look in which temperature region the particle surfaces and grid box region
    # released
    for i in range(le):
        
        lo = find_down(Lons,lons0[i]); lon = Lons[lo];
        la = find_down(Lats,lats0[i]); lat = Lats[la];
        
        j = np.where(np.logical_and(vLons==lon,vLats==lat))[0][0]  
           
        tssalin[j] = np.append(tssalin[j], salin[i])
        tstemp[j] = np.append(tstemp[j], temp[i])   
        tslon[j] = np.append(tslon[j], lonadv[i])
        tslat[j] = np.append(tslat[j], latadv[i])  
        tslon0[j] = np.append(tslon0[j], lons0[i])
        tslat0[j] = np.append(tslat0[j], lats0[i])  
        tsage[j] = np.append(tsage[j], ageadv[i])   


    for j in range(len(vLons)):
        maxlents = max(maxlents,len(tstemp[j]))                    
                    
    return tstemp, tssalin, tslon, tslat, tslon0, tslat0, tsage, maxlents

def construct_tsfix(tsfixtemp, tsfixsalin, fixlon, fixlat,  vLons, vLats, Lons, Lats, fixtemp, fixsalin):
    """Return a list with first dimension len(VLons) for all the bottom grid boxes
    and second dimension the temperatures at the surface related to the bottom box"""  
    sle = len(fixlon)
    
    maxlentsfix = 0
    
    for i in range(sle):
        lo = find_down(Lons,fixlon[i]); lon = Lons[lo];
        la = find_down(Lats,fixlat[i]); lat = Lats[la];
        
        j = np.where(np.logical_and(vLons==lon,vLats==lat))[0][0]

        tsfixtemp[j] = np.append(tsfixtemp[j],fixtemp[i,:])     
        tsfixsalin[j] = np.append(tsfixsalin[j],fixsalin[i,:])
        
    for i in range(len(vLons)):
        maxlentsfix = max(maxlentsfix,len(tsfixtemp[i]))  
                    
    return tsfixtemp, tsfixsalin, maxlentsfix

lat0 = np.array(pfile['lat0'][:])
lon0 = np.array(pfile['lon0'][:])
latadv = np.array(pfile['lat'][:])
lonadv = np.array(pfile['lon'][:])
ageadv = np.array(pfile['age'][:])
lonfix = pfilefix['lon'][:,0]
latfix = pfilefix['lat'][:,0]
temp = np.array(pfile['temp'][:])
tempfix = pfilefix['temp'][:,:]
salin = np.array(pfile['salin'][:])
salinfix = pfilefix['salin'][:,:]

print 'full fixlat: ', latfix

print 'any nans:   ',np.isnan(lat0).any(), np.isnan(lon0).any(), np.isnan(lonfix).any(), np.isnan(latfix).any(), np.isnan(temp).any(), np.isnan(tempfix).any()
print 'all shapes: ', lat0.shape, lon0.shape, temp.shape, lonfix.shape, latfix.shape, tempfix.shape
print vLons.shape, vLats.shape, Lons.shape, Lats.shape

start = time.time() 
tstemp = np.empty((len(vLons),), dtype=object)  
for i in range(len(tstemp)):tstemp[i] = [];
tssalin = np.empty((len(vLons),), dtype=object)  
for i in range(len(tssalin)):tssalin[i] = [];
tslon = np.empty((len(vLons),), dtype=object)  
for i in range(len(tslon)):tslon[i] = [];
tslat = np.empty((len(vLons),), dtype=object)  
for i in range(len(tslat)):tslat[i] = [];
tslon0 = np.empty((len(vLons),), dtype=object)  
for i in range(len(tslon0)):tslon0[i] = [];
tslat0 = np.empty((len(vLons),), dtype=object)  
for i in range(len(tslat0)):tslat0[i] = [];
tsage = np.empty((len(vLons),), dtype=object)  
for i in range(len(tsage)):tsage[i] = [];
    
tstemp, tssalin, tslon, tslat, tslon0, tslat0, tsage, maxlents =  construct_ts(tstemp, tssalin, tslon, tslat, tslon0, tslat0, tsage, lat0, lon0, vLons, vLats, Lons, Lats, temp, salin, lonadv, latadv, ageadv)
print 'time \'ts\' (minutes): ', ((time.time()-start)/60.)

ts = np.array(ts)

start = time.time()   
tsfixtemp = np.empty((len(vLons),), dtype=object)  
for i in range(len(tsfixtemp)):tsfixtemp[i] = [];
tsfixsalin = np.empty((len(vLons),), dtype=object)  
for i in range(len(tsfixsalin)):tsfixsalin[i] = [];
    
tsfixtemp, tsfixsalin, maxlentsfix =  construct_tsfix(tsfixtemp, tsfixsalin, lonfix,latfix,  vLons, vLats, Lons, Lats, tempfix, salinfix)

print 'time \'tsfix\' (minutes): ', ((time.time()-start)/60.)     
print 'maxlents:', maxlents
print 'maxlentsfix:', maxlentsfix
  
if(type(sp)==str):
    if(sp=='s1'):
        dataset = Dataset(dirWrite  + 'timeseries_per_location_inclatlonadv_ddeg%d_spincreaseS1_dd%d.nc'%(ddeg,int(dd)),'w',format='NETCDF4_CLASSIC')
    else:
        dataset = Dataset(dirWrite  + 'timeseries_per_location_inclatlonadv_ddeg%d_spincreaseS2_dd%d.nc'%(ddeg,int(dd)),'w',format='NETCDF4_CLASSIC')
else:
    dataset = Dataset(dirWrite  + 'timeseries_per_location_inclatlonadv_ddeg%d_sp%d_dd%d.nc'%(ddeg, int(sp),int(dd)),'w',format='NETCDF4_CLASSIC')

traj = dataset.createDimension('tslen', maxlents)
trajfix = dataset.createDimension('tslenfix', maxlentsfix)
vlatlon = dataset.createDimension('vLons', len(vLons))
lons = dataset.createDimension('Lons', len(Lons))
lats = dataset.createDimension('Lats', len(Lats))

vLon = dataset.createVariable('vLons', np.float64, ('vLons',))
vLat = dataset.createVariable('vLats', np.float64, ('vLons',))
Lon = dataset.createVariable('Lons', np.float64, ('Lons',))
Lat = dataset.createVariable('Lats', np.float64, ('Lats',))

temps = dataset.createVariable('temp', np.float64, ('vLons','tslen',))
salins = dataset.createVariable('salin', np.float64, ('vLons','tslen',))
lat = dataset.createVariable('lat', np.float64, ('vLons','tslen',))
lon = dataset.createVariable('lon', np.float64, ('vLons','tslen',))
lat0 = dataset.createVariable('lat0', np.float64, ('vLons','tslen',))
lon0 = dataset.createVariable('lon0', np.float64, ('vLons','tslen',))
ages = dataset.createVariable('age', np.float64, ('vLons','tslen',))

fixtemps = dataset.createVariable('fixtemp', np.float64, ('vLons','tslenfix',))
fixsalins = dataset.createVariable('fixsalin', np.float64, ('vLons','tslenfix',))

vLon[:] = vLons
vLat[:] = vLats
Lon[:] = Lons
Lat[:] = Lats

for j in range(len(vLons)):
    le = len(tstemp[j])
    if(le>0):
        temps[j,:le] = np.array(tstemp[j])
        salins[j,:le] = np.array(tssalin[j])
        lon[j,:le] = np.array(tslon[j])
        lat[j,:le] = np.array(tslat[j])
        lon0[j,:le] = np.array(tslon0[j])
        lat0[j,:le] = np.array(tslat0[j])
        ages[j,:le] = np.array(tsage[j])

    fle = len(tsfixtemp[j])
    if(fle>0):
        fixtemps[j,:fle] = np.array(tsfixtemp[j])
        fixsalins[j,:fle] = np.array(tsfixsalin[j])

dataset.close()
