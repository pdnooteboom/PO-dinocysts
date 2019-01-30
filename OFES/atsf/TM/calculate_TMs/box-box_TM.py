# -*- coding: utf-8 -*-
"""
Created on Tue May 29 10:03:01 2018

@author: nooteboom
"""
import numpy as np
from netCDF4 import Dataset
import Create_TM as tm

dd = 10. #Dwelling depth
sp = 11. #Sinking speed
res = 1 #Degree resolution on which the particles are released

dirRead = '/projects/0/palaeo-parcels/OFESres/particledata/'
dirWrite = '/projects/0/palaeo-parcels/OFESres/TM/box-box/'

ddeg = 1 # Binning size
Lons = np.arange(0,360,ddeg)
Lats = np.arange(-75,75,ddeg)

lats = np.array([])
lons = np.array([])
lats0 = np.array([])
lons0 = np.array([])
for posidx in range(108):

    if(posidx==84 or posidx==85):
        print'no posidx file:'
        print posidx
    else:
        print 'posidx: ',posidx
        pfile = Dataset(dirRead + 'sp%d_dd%d/grid_id'%(int(sp),int(dd)) +str(posidx)+'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res) + ".nc")
        #pfile_surf = Dataset(dirRead + 'surface/surfacegrid_id'+str(posidx)+'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res) )    
        
        lats = np.append(lats,pfile['lat'][:])
        lats0 = np.append(lats0,pfile['lat0'][:])
        lons = np.append(lons,pfile['lon'][:])
        if(pfile['lon0'][0]<0):
            lons0 = np.append(lons0,pfile['lon0'][:]+360)   
        else:
            lons0 = np.append(lons0,pfile['lon0'][:])  
        print pfile['lon0'][0]
        print pfile['lon'][0]


    
lon=np.concatenate((np.expand_dims(lons0[:],axis=1),np.expand_dims(lons[:],axis=1)),axis=1)
lat=np.concatenate((np.expand_dims(lats0[:],axis=1),np.expand_dims(lats[:],axis=1)),axis=1)   
    
P = tm.TM_ParticleSet(lon,lat, Lons, Lats)
P.setup_TM(ddeg,Lons,Lats)    
    
    
P.save_TM(dirWrite+'TMglobal_bin'+str(ddeg)+'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res) )
    