# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 13:48:07 2019

@author: nooteboom
"""

import numpy as np
import math


def find_down(array,value):
    if(value>=array[0]):
        array = [n-value for n in array]
        idx = np.array([n for n in array if n<0]).argmax()
    else:
        idx = np.nan
    return idx 

def region_index(lats0, lons0):
    """Return the indices of the grid box at the bottom (Bi)""" 
    Bi = np.zeros(lons0.shape[0]) #The grid cell index where a particle is released at the bottom

    for i in range(lons0.shape[0]):
        lo = find_down(Lons,lons0[i]);lon = Lons[lo];
        la = find_down(Lats,lats0[i]);lat = Lats[la];
        Bi[i] = np.where(np.logical_and(vLons==lon,vLats==lat))[0]                  
    return Bi
#%%
ddeg = 2 # resolution of the binning
spl = [3, 6, 's1','s2', 11, 25, 50, 100, 200, 500]
dd = 10
res = 1

tmdir = '/Users/nooteboom/Documents/PhD/parcels/OFES_global/Transtition_matrices/'

TM = []
for sp in spl:
    if(type(sp)==str):
        data = np.load(tmdir + 'output/box-box/TMglobal_bin'+str(ddeg)+'_dd'+str(int(dd)) +'_sp'+sp+"_res"+str(res) + '.npz')
    else:
        data = np.load(tmdir + 'output/box-box/TMglobal_bin'+str(ddeg)+'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res) + '.npz')
    TM.append(data['TM'][:])

Lons = data['Lons'][:]
Lats = data['Lats'][:]
Lonss, Latss = np.meshgrid(Lons, Lats)

vLons = Lonss.reshape(Lonss.size)
vLats = Latss.reshape(Latss.size)

maxlat = 70
latind = np.where(vLats<maxlat)

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

vLons, vLats = np.meshgrid(Lons, Lats)
vLons = vLons.flatten(); vLats = vLats.flatten();

surface = np.full(vLons.shape[0], 0)

for i in range(len(vLons)):
    lonmin = vLons[i]; lonmax = vLons[i]+ddeg;
    latmin = vLats[i]; latmax = vLats[i]+ddeg;   
    dis1 = distance((latmin,lonmin),(latmax,lonmin))
    dis2 = distance((latmin,lonmin),(latmin, lonmax))
    surface[i] = dis1 * dis2
    
#%% boxplot of the zonneveld locations
zonread = '/Users/nooteboom/Documents/PhD/parcels/secondpart/'

zondata = np.load(zonread + 'Zonneveldt2013_data.npy').item()
zonlat = zondata['data'][:,1].astype(float)
zonlon = zondata['data'][:,2].astype(float)

zonlon[zonlon<0] = zonlon[zonlon<0] + 360
zonlon = zonlon[zonlat<min(maxlat,75)]
zonlat = zonlat[zonlat<min(maxlat,75)]

zon_in_TM = region_index(zonlat, zonlon).astype(np.int)

zondata = []
for i in range(len(TM)):
    connected = (TM[i]>0).astype(int)
    for j in range(connected.shape[0]):
        connected[j] = connected[j] * surface
    zondata.append(np.sum(connected[:,(zon_in_TM)], axis=0)/1000000.)


zondata_diag = [np.diagonal(TM[i])[(zon_in_TM)][np.diagonal(TM[i])[(zon_in_TM)]>0]*100 for i in range(len(spl))]

#%% boxplot of the data: 
plotdata = []
for i in range(len(TM)):
    connected = (TM[i][:,latind[0]]>0).astype(int)
    for j in range(connected.shape[0]):
        connected[j] = connected[j] * surface[latind]
    plotdata.append(np.sum(connected, axis=0)/1000000.)#/surface.astype(float))    

plotdata_diag = [np.diagonal(TM[i])[np.diagonal(TM[i])>0]*100 for i in range(len(spl))]

areas = np.array([])
diagonal = np.array([])
sinkingspeed = np.array([])
locations = np.array([])
sinkingspeed_diag = np.array([])
locations_diag = np.array([])
for i in range(len(spl)):
    
    if(spl[i]=='s1'):
        name = 'SC1'
    elif(spl[i]=='s2'):
        name = 'SC2'
    else:
        name = str(spl[i])    
    
    diagonal = np.append(diagonal,plotdata_diag[i]); diagonal = np.append(diagonal,zondata_diag[i]); 
    sinkingspeed_diag = np.append(sinkingspeed_diag,np.full(len(plotdata_diag[i]),name)); sinkingspeed_diag = np.append(sinkingspeed_diag,np.full(len(zondata_diag[i]),name)); 
    locations_diag = np.append(locations_diag,np.full(len(plotdata_diag[i]),'global')); locations_diag = np.append(locations_diag,np.full(len(zondata_diag[i]),'core top sample sites'));     


    areas = np.append(areas,plotdata[i]); areas = np.append(areas,zondata[i]); 
    sinkingspeed = np.append(sinkingspeed,np.full(len(plotdata[i]),name)); sinkingspeed = np.append(sinkingspeed,np.full(len(zondata[i]),name)); 
    locations = np.append(locations,np.full(len(plotdata[i]),'global')); locations = np.append(locations,np.full(len(zondata[i]),'core top sample sites'));     

#%% The average distance part

#%Load avg drift distance
datadist = np.load(tmdir + 'output/box-avgdist/TM_box-avgdist_ddeg%d_sp'%(ddeg)+str(sp)+'_dd%d.npz'%( dd))
TMd = datadist['TM'][:]
Lonsd = datadist['Lons'][:]
Latsd = datadist['Lats'][:]
vLons, vLats = np.meshgrid(Lonsd, Latsd)
vLons = vLons.flatten(); vLats = vLats.flatten();

Lonsd[Lonsd>=360] -= 360; 

TMd = TMd[np.logical_and(vLons<360,vLats<=73)]

TMd = []
for sp in spl:
    if(type(sp)==str):
        data = np.load(tmdir + 'output/box-avgdist/TM_box-avgdist_ddeg%d_sp'%(ddeg)+sp+'_dd%d.npz'%(dd))
    else:    
        data = np.load(tmdir + 'output/box-avgdist/TM_box-avgdist_ddeg%d_sp%d_dd%d.npz'%(ddeg, sp, dd))
    TMd.append(data['TM'][:])
    
distance = np.array([])
sinkingspeed_dist = np.array([])
locations_dist = np.array([])
for i in range(len(spl)):
    distance = np.append(distance,TMd[i]); distance = np.append(distance,TMd[i][(zon_in_TM)]); 
    if(spl[i]=='s1'):
        name = 'SC1'
    elif(spl[i]=='s2'):
        name = 'SC2'
    else:
        name = str(spl[i])
    sinkingspeed_dist = np.append(sinkingspeed_dist,np.full(len(TMd[i]),name)); sinkingspeed_dist = np.append(sinkingspeed_dist,np.full(len(TMd[i][(zon_in_TM)]),name)); 
    locations_dist = np.append(locations_dist,np.full(len(TMd[i]),'global')); locations_dist = np.append(locations_dist,np.full(len(TMd[i][(zon_in_TM)]),'core top sample sites'));     
 

np.savez('toplot_TM-boxplot', diagonal = diagonal, sinkingspeed_diag = sinkingspeed_diag, 
         locations_diag=locations_diag, areas=areas, sinkingspeed=sinkingspeed, locations=locations,
         distance = distance, sinkingspeed_dist=sinkingspeed_dist, locations_dist = locations_dist)






