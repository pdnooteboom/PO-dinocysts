# -*- coding: utf-8 -*-
"""
Created on Tue May 29 10:03:01 2018

@author: nooteboom
"""
import numpy as np
from netCDF4 import Dataset
import Create_TM as tm

dd = 10. #Dwelling depth
sp = 6. #Sinking speed
res = 1 #Degree resolution on which the particles are released

OFEScompare = False

dirRead = '/NEMO/NEMOres/atsf/particlefiles/'
dirWrite = '/NEMO/NEMOres/atsf/TM/'

ddeg = 2 # Binning size in degrees
if(not OFEScompare):
    Lons = np.arange(-180,180,ddeg)
    Lats = np.arange(-77,85,ddeg)
else:
    Lons = np.arange(0,360,ddeg)
    Lats = np.arange(-75,75,ddeg)

lats = np.array([])
lons = np.array([])
lats0 = np.array([])
lons0 = np.array([])
for posidx in range(17):    
    
    pfile = Dataset(dirRead + 'global_sp%d_dd%d/grid_nobgc_id'%(int(sp),int(dd))+str(posidx)+'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res)+'.nc')  
    
    lats = np.append(lats,pfile['lat'][:])
    lats0 = np.append(lats0,pfile['lat0'][:])
    lons = np.append(lons,pfile['lon'][:])
    lons0 = np.append(lons0,pfile['lon0'][:])     
    
nanv = np.where(np.logical_and(~np.isnan(lats),~np.isnan(lons)))
lats = lats[nanv]
lons = lons[nanv]
lats0 = lats0[nanv]
lons0 = lons0[nanv]
    
lon=np.concatenate((np.expand_dims(lons0[:],axis=1),np.expand_dims(lons[:],axis=1)),axis=1)
lat=np.concatenate((np.expand_dims(lats0[:],axis=1),np.expand_dims(lats[:],axis=1)),axis=1)   

if(OFEScompare):
    lon[np.where(lon<0)] += 360
    
P = tm.TM_ParticleSet(lon,lat,ddeg, Lons, Lats)
P.setup_TM(ddeg,Lons,Lats)

if(OFEScompare):
    P.save_TM(dirWrite+'box-box/TMglobal_compareOFES_bin'+str(ddeg)+'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res) )
else:    
    P.save_TM(dirWrite+'box-box/TMglobal_bin'+str(ddeg)+'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res) )
    
