# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 13:16:25 2018

@author: nooteboom
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 16:41:17 2018

@author: nooteboom
"""

import numpy as np
import matplotlib.pylab as plt
import matplotlib
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Polygon
import seaborn as sns
import pandas as pd

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
font = {'family' : 'Helvetica',
#'weight' : 'bold',
'size'   : 60}

matplotlib.rc('font', **font)

# True if surface instead of amount of surface grid boxes
surfacebool = True

ddeg = 2 # resolution of the binning
spl = [6, 11, 50]
dd = 10
res = 1

#Set for Pdensity plots:
scale = 'exponential'
k_depth = 'proportion'
outlier_prop = 0.5
whis = [5, 95]

meanpointprops = dict(marker='D', markeredgecolor='black',
                      markerfacecolor='firebrick', markersize=7)

tmdir = '/Users/nooteboom/Documents/PhD/parcels/NEMO/atsf/Transition_matrices/'

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
import math

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

#%% boxplot of the zonneveld locations
zonread = '/Users/nooteboom/Documents/PhD/parcels/secondpart/'

zondata = np.load(zonread + 'Zonneveldt2013_data.npy').item()
zonlat = zondata['data'][:,1].astype(float)
zonlon = zondata['data'][:,2].astype(float)

zonlon[zonlon<0] = zonlon[zonlon<0] + 360
zonlon = zonlon[zonlat<min(maxlat,75)]
zonlat = zonlat[zonlat<min(maxlat,75)]

zon_in_TM = region_index(zonlat, zonlon).astype(np.int)

if(surfacebool):
    zondata = []
    for i in range(len(TM)):
        connected = (TM[i]>0).astype(int)
        for j in range(connected.shape[0]):
            connected[j] = connected[j] * surface
        zondata.append(np.sum(connected[:,(zon_in_TM)], axis=0)/1000000.)#/surface[(zon_in_TM)].astype(float))
else:
    zondata = [np.sum(TM[i][:,(zon_in_TM)]>0, axis=0) for i in range(len(spl))]

zondata_diag = [np.diagonal(TM[i])[(zon_in_TM)][np.diagonal(TM[i])[(zon_in_TM)]>0]*100 for i in range(len(spl))]
#zondata_diag = [np.diagonal(TM[i])[(zon_in_TM)]*100 for i in range(len(spl))]
#%% boxplot of the data: 

if(surfacebool):
    plotdata = []
    for i in range(len(TM)):
        connected = (TM[i][:,latind[0]]>0).astype(int)
        for j in range(connected.shape[0]):
            connected[j] = connected[j] * surface[latind]
        plotdata.append(np.sum(connected, axis=0)/1000000.)#/surface.astype(float))    
else:
    plotdata = [np.sum(TM[i]>0, axis=0) for i in range(len(spl))]#, np.sum(TM[3]>0, axis=0)]

plotdata_diag = [np.diagonal(TM[i])[np.diagonal(TM[i])>0]*100 for i in range(len(spl))]
#plotdata_diag = [np.diagonal(TM[i])*100 for i in range(len(spl))]

#or i in range(len(spl)):
#    print 'perc: ',np.percentile(np.diagonal(TM[i])*100,20),'    sp:',spl[i]
#%% boxplot of different sinking velocities


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
    locations_diag = np.append(locations_diag,np.full(len(plotdata_diag[i]),'all')); locations_diag = np.append(locations_diag,np.full(len(zondata_diag[i]),'measured cores'));     


    areas = np.append(areas,plotdata[i]); areas = np.append(areas,zondata[i]); 
    sinkingspeed = np.append(sinkingspeed,np.full(len(plotdata[i]),name)); sinkingspeed = np.append(sinkingspeed,np.full(len(zondata[i]),name)); 
    locations = np.append(locations,np.full(len(plotdata[i]),'all')); locations = np.append(locations,np.full(len(zondata[i]),'measured cores'));     

if(surfacebool):
    d = {'$10^6\ km^2$':areas,
     'sinking speed ($m\ day^{-1}$)':sinkingspeed,
    'locations':locations }    
else:
    d = {'# boxes':areas,
     'sinking speed ($m\ day^{-1}$)': sinkingspeed,
    'locations': locations }
                                                                         
np.save('TM_boxplot_secondsubplot',d) 


df = pd.DataFrame(data=d)
sns.set(style="whitegrid", font='Helvetica', font_scale=1.2)

fig = plt.figure(figsize=(15,10))

#fig.subplots_adjust(bottom=0.25)

#plt.subplot(124)
ax2 = plt.subplot2grid((2,19),(0,10), colspan=9)


plt.title('(b)', fontsize=25)
order = [str(i) for i in spl];
for i in range(len(order)): 
    if(order[i]=='s1'): 
        order[i] = 'SC1'
    elif(order[i]=='s2'):
        order[i] = 'SC2'
order = np.array(order)
if(surfacebool):
#    ax = sns.boxplot(x='sinking speed ($m\ day^{-1}$)', y='area surface ($10^6\ km^2$)', hue='locations', data=df, order=order)#, showfliers=False)    
#    ax = sns.boxenplot(x='sinking speed ($m\ day^{-1}$)', y='area surface ($10^6\ km^2$)', hue='locations', data=df, order=order, k_depth=k_depth, outlier_prop=outlier_prop, scale=scale)
    ax = sns.boxplot(x='sinking speed ($m\ day^{-1}$)', y='$10^6\ km^2$', hue='locations', data=df, showmeans=True, meanprops=meanpointprops, order=order, whis=whis,  showfliers=True)
  #  ax = sns.violinplot(x='sinking speed ($m\ day^{-1}$)', y='area surface ($10^6\ km^2$)', hue='locations', data=df, order=order)
else:
    ax = sns.boxplot(x='sinking speed ($m\ day^{-1}$)', y='# boxes', hue='locations', data=df, order=order)#, showfliers=False)
#ax.set(yscale= 'log')


#plt.show()


#ax.set(yscale= 'log')

#
#ax = sns.stripplot(x='sinking speed (m/day)', y='# boxes', hue='locations', data=df)
##ax.set(yscale= 'log')
#plt.show()

#ax = sns.swarmplot(x='sinking speed (m/day)', y='# boxes', hue='locations', data=df, dodge=True)
#ax.set(yscale= 'log')
#plt.show()


#% Boxplot of the diagonal of the transition matrix


d =     {'%':diagonal,
         'sinking speed ($m\ day^{-1}$)': sinkingspeed_diag,
        'locations': locations_diag}
        
np.save('TM_boxplot_firstsubplot',d)

df = pd.DataFrame(data=d)

sns.set(style="whitegrid", font='Helvetica', font_scale=1.2)

#plt.subplot(121)
ax2 = plt.subplot2grid((2,19),(0,0), colspan=9)
#ax = sns.boxenplot(x='sinking speed (m/day)', y='(%)', hue='locations', data=df, order=order, k_depth=k_depth, outlier_prop=outlier_prop, scale=scale)#,  showfliers=True)
ax = sns.boxplot(x='sinking speed ($m\ day^{-1}$)', y='%', hue='locations', data=df, showmeans=True, meanprops=meanpointprops, order=order, whis=whis,  showfliers=False)
ax.legend_.remove()
#ax.set(yscale= 'log')
plt.title('(a)', fontsize=25)

#plt.show()

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
    
diagonal = np.array([])
sinkingspeed_diag = np.array([])
locations_diag = np.array([])
for i in range(len(spl)):
    diagonal = np.append(diagonal,TMd[i]); diagonal = np.append(diagonal,TMd[i][(zon_in_TM)]); 
    if(spl[i]=='s1'):
        name = 'SC1'
    elif(spl[i]=='s2'):
        name = 'SC2'
    else:
        name = str(spl[i])
    sinkingspeed_diag = np.append(sinkingspeed_diag,np.full(len(TMd[i]),name)); sinkingspeed_diag = np.append(sinkingspeed_diag,np.full(len(TMd[i][(zon_in_TM)]),name)); 
    locations_diag = np.append(locations_diag,np.full(len(TMd[i]),'all')); locations_diag = np.append(locations_diag,np.full(len(TMd[i][(zon_in_TM)]),'measured cores'));     


d =     {'km':diagonal,
         'sinking speed ($m\ day^{-1}$)': sinkingspeed_diag,
        'locations': locations_diag}

df = pd.DataFrame(data=d)        
ax2 = plt.subplot2grid((2,19),(1,5), colspan=9)
#ax = sns.boxenplot(x='sinking speed (m/day)', y='(%)', hue='locations', data=df, order=order, k_depth=k_depth, outlier_prop=outlier_prop, scale=scale)#,  showfliers=True)
ax = sns.boxplot(x='sinking speed ($m\ day^{-1}$)', y='km', hue='locations', data=df, showmeans=True, meanprops=meanpointprops, meanline=False, order=order, whis=whis,  showfliers=False)
ax.legend_.remove()
#ax.set(yscale= 'log')
plt.title('(c)', fontsize=25)       
plt.savefig('/Users/nooteboom/Documents/PhD/firstpaper/articleplots/nemoplots/' + 'boxplot_TMsinksensitivity_withavgdist.pdf', bbox_inches="tight")
plt.show()
