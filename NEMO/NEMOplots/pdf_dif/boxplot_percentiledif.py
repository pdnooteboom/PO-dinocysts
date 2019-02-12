# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 13:16:25 2018

@author: nooteboom
"""

import numpy as np
import matplotlib.pylab as plt
import matplotlib
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
'size'   : 40}

matplotlib.rc('font', **font)

ddeg = 1 # resolution of the binning
spl = [6,11,50]
dd = 10
res = 1

#Set for Pdensity plots:
scale = 'linear'
k_depth = 'proportion'
outlier_prop = 0.2
whis = [1,99]#[5,95]#

order = [str(i) for i in spl];
for i in range(len(order)): 
    if(order[i]=='s1'): 
        order[i] = 'SC1'
    elif(order[i]=='s2'):
        order[i] = 'SC2'
order = np.array(order)

tmdir = '/Users/nooteboom/Documents/PhD/parcels/NEMO/atsf/Transition_matrices/'

var = 'temp'
var2 = 'salin'

n = 5

low_var1 = []
high_var1 = []
low_var2 = []
high_var2 = []
for sp in spl:
    if(type(sp)==str):
        datalow = np.load(tmdir + 'output/box-percentile/TM_box-%dpercentiledif_'%(n)+var+'_ddeg%d_sp'%(ddeg)+sp+'_dd%d'%(int(dd))+'.npz')
        low_var1.append(datalow['TM'][:])
        datahigh = np.load(tmdir + 'output/box-percentile/TM_box-%dpercentiledif_'%(100-n)+var+'_ddeg%d_sp'%(ddeg)+sp+'_dd%d'%(int(dd))+'.npz')
        high_var1.append(datahigh['TM'][:])
    
        datalow = np.load(tmdir + 'output/box-percentile/TM_box-%dpercentiledif_'%(n)+var2+'_ddeg%d_sp'%(ddeg)+sp+'_dd%d'%(int(dd))+'.npz')
        low_var2.append(datalow['TM'][:])
        datahigh = np.load(tmdir + 'output/box-percentile/TM_box-%dpercentiledif_'%(100-n)+var2+'_ddeg%d_sp'%(ddeg)+sp+'_dd%d'%(int(dd))+'.npz')
        high_var2.append(datahigh['TM'][:])        
    else:
        datalow = np.load(tmdir + 'output/box-percentile/TM_box-%dpercentiledif_'%(n)+var+'_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd))+'.npz')
        low_var1.append(datalow['TM'][:])
        datahigh = np.load(tmdir + 'output/box-percentile/TM_box-%dpercentiledif_'%(100-n)+var+'_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd))+'.npz')
        high_var1.append(datahigh['TM'][:])
    
        datalow = np.load(tmdir + 'output/box-percentile/TM_box-%dpercentiledif_'%(n)+var2+'_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd))+'.npz')
        low_var2.append(datalow['TM'][:])
        datahigh = np.load(tmdir + 'output/box-percentile/TM_box-%dpercentiledif_'%(100-n)+var2+'_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd))+'.npz')
        high_var2.append(datahigh['TM'][:])


Lons = datalow['Lons'][:]
Lats = datalow['Lats'][:]
Lonss, Latss = np.meshgrid(Lons, Lats)

vLons = Lonss.reshape(Lonss.size)
vLats = Latss.reshape(Latss.size)

maxlat = 70

#%% boxplot of the zonneveld locations
zonread = '/Users/nooteboom/Documents/PhD/parcels/secondpart/'

zondata = np.load(zonread + 'Zonneveldt2013_data.npy').item()
zonlat = zondata['data'][:,1].astype(float)
zonlon = zondata['data'][:,2].astype(float)

zonlon[zonlon<0] = zonlon[zonlon<0] + 360
zonlon = zonlon[zonlat<min(maxlat,75)]
zonlat = zonlat[zonlat<min(maxlat,75)]

zon_in_TM = region_index(zonlat, zonlon).astype(np.int)

zondata_low_var1 = [low_var1[i][(zon_in_TM)] for i in range(len(spl))]
zondata_low_var2 = [low_var2[i][(zon_in_TM)] for i in range(len(spl))]
zondata_high_var1 = [high_var1[i][(zon_in_TM)] for i in range(len(spl))]
zondata_high_var2 = [high_var2[i][(zon_in_TM)] for i in range(len(spl))]
#%% boxplot of the data: 


plotdata_low_var1 = [low_var1[i] for i in range(len(spl))]
plotdata_low_var2 = [low_var2[i] for i in range(len(spl))]
plotdata_high_var1 = [high_var1[i] for i in range(len(spl))]
plotdata_high_var2 = [high_var2[i] for i in range(len(spl))]
#%% boxplot of different sinking velocities


d_low_var2 = np.array([])
d_low_var1 = np.array([])
d_high_var1 = np.array([])
d_high_var2 = np.array([])
sinkingspeed = np.array([])
locations = np.array([])
for i in range(len(spl)):
    d_low_var1 = np.append(d_low_var1,plotdata_low_var1[i]); d_low_var1 = np.append(d_low_var1,zondata_low_var1[i]); 
    d_low_var2 = np.append(d_low_var2,plotdata_low_var2[i]); d_low_var2 = np.append(d_low_var2,zondata_low_var2[i]); 
    d_high_var1 = np.append(d_high_var1,plotdata_high_var1[i]); d_high_var1 = np.append(d_high_var1,zondata_high_var1[i]); 
    d_high_var2 = np.append(d_high_var2,plotdata_high_var2[i]); d_high_var2 = np.append(d_high_var2,zondata_high_var2[i]);     

    if(spl[i]=='s1'):
        name = 'SC1'
    elif(spl[i]=='s2'):
        name = 'SC2'
    else:
        name = str(spl[i])
    
    sinkingspeed = np.append(sinkingspeed,np.full(len(plotdata_low_var1[i]),name)); sinkingspeed = np.append(sinkingspeed,np.full(len(zondata_low_var1[i]),name)); 
    locations = np.append(locations,np.full(len(plotdata_low_var1[i]),'all')); locations = np.append(locations,np.full(len(zondata_low_var1[i]),'measured cores'));  
    
#%%
plt.figure(figsize=(17,8))
    
#%% Boxplot 5% percentile var1


d =     {'$^{\circ}$C':d_low_var1,
         'sinking speed ($m\ day^{-1}$)': sinkingspeed,
        'locations': locations}
df = pd.DataFrame(data=d)

sns.set(style="whitegrid", font='Helvetica', font_scale=1.2)

plt.subplot(221)
#ax = sns.boxenplot(x='sinking speed (m/day)', y='$^{\circ}$C', hue='locations', data=df, order=order, k_depth=k_depth, scale=scale,outlier_prop=outlier_prop)#, linewidth=2), notch=True,  showfliers=False, whis=3
ax = sns.boxplot(x='sinking speed ($m\ day^{-1}$)', y='$^{\circ}$C', hue='locations', data=df, order=order, whis=whis,  showfliers=False)
#ax = sns.violinplot(x='sinking speed (m/day)', y='$^{\circ}$C', hue='locations', data=df, order=order)
#plt.ylim(-30,18)
ax.set_xlabel('')
ax.set(xticklabels=[])
#ax.set(yscale= 'log')
plt.title('(a)', fontsize=25)    

#%% Boxplot 95% percentile var1

d = {'$^{\circ}$C':d_high_var1,
 'sinking speed ($m\ day^{-1}$)': sinkingspeed,
'locations': locations }
    
df = pd.DataFrame(data=d)
sns.set(style="whitegrid", font='Helvetica', font_scale=1.2)

plt.subplot(222)
plt.title('(b)', fontsize=25)

#ax = sns.boxenplot(x='sinking speed (m/day)', y='$^{\circ}$C', hue='locations', data=df, order=order, k_depth=k_depth, scale=scale,outlier_prop=outlier_prop)
ax = sns.boxplot(x='sinking speed ($m\ day^{-1}$)', y='$^{\circ}$C', hue='locations', data=df, order=order, whis=whis,  showfliers=False)
#plt.ylim(-30,18)
ax.set(xticklabels=[])
ax.legend_.remove()
ax.set_ylabel('')
ax.set_xlabel('')
#%% Boxplot 5% percentile var2


d =     {'$^{\circ}$C':d_low_var2,
         'sinking speed ($m\ day^{-1}$)': sinkingspeed,
        'locations': locations}
df = pd.DataFrame(data=d)

sns.set(style="whitegrid", font='Helvetica', font_scale=1.2)

plt.subplot(223)
#ax = sns.boxenplot(x='sinking speed (m/day)', saturation=1, y='$^{\circ}$C', hue='locations', data=df, order=order, k_depth=k_depth, scale=scale,outlier_prop=outlier_prop)#, linewidth=2), notch=True,  showfliers=False, whis=3
ax = sns.boxplot(x='sinking speed ($m\ day^{-1}$)', y='$^{\circ}$C', hue='locations', data=df, order=order, whis=whis,  showfliers=False)
ax.set_ylabel('PSU')
#ax.set(yscale= 'log')
plt.title('(c)', fontsize=25)    
#plt.ylim(-30,7)
ax.legend_.remove()

#%% Boxplot 95% percentile var2

d = {'$^{\circ}$C':d_high_var2,
 'sinking speed ($m\ day^{-1}$)': sinkingspeed,
'locations': locations }
    
df = pd.DataFrame(data=d)
sns.set(style="whitegrid", font='Helvetica', font_scale=1.2)

plt.subplot(224)
plt.title('(d)', fontsize=25)

#ax = sns.boxenplot(x='sinking speed (m/day)', y='$^{\circ}$C', hue='locations', data=df, order=order, k_depth=k_depth, scale=scale,outlier_prop=outlier_prop)
ax = sns.boxplot(x='sinking speed ($m\ day^{-1}$)', y='$^{\circ}$C', hue='locations', data=df, order=order, whis=whis,  showfliers=False)
#plt.ylim(-2,6)
ax.legend_.remove()
ax.set_ylabel('')

#%%
plt.savefig('/Users/nooteboom/Documents/PhD/firstpaper/articleplots/nemoplots/' + 'boxplotNEMO_PDFsinksensitivity.pdf', bbox_inches="tight")
plt.show()

