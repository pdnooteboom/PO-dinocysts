# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 09:56:50 2019

@author: nooteboom
"""

import numpy as np 
from sklearn.cluster import KMeans
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import seaborn as sns
from copy import copy
import matplotlib

#%% Load Zonneveld data:
zonread = '/Users/nooteboom/Documents/PhD/parcels/secondpart/'

zondata = np.load(zonread + 'Zonneveldt2013_data.npy').item()

species = zondata['variables'][19:-2]

maxlat = -12

sdata = zondata['data'][:,19:-2][zondata['data'][:,1].astype(float)<maxlat].astype(float)

lats = zondata['data'][:,1][zondata['data'][:,1].astype(float)<maxlat].astype(float)
lons = zondata['data'][:,2][zondata['data'][:,1].astype(float)<maxlat].astype(float)
lons[lons<0] += 360

#%% k-means clustering
n_clusters = 7

kmeans = KMeans(n_clusters=n_clusters, random_state=2).fit(sdata)
czon = copy(kmeans.labels_)
colors = ['navy', 'royalblue', 'lightskyblue', 'seagreen','yellow','orange','red','darkred']
symbols = ['D','s','v','o','d','*','+','p']
#%%Load the model PDFs after advection
var = 'temp'
ddeg = 1
sp = 6
dd = 10

tail_perc = 5

dirRead = '/OFES/OFESres/TM/input/'

pf = Dataset(dirRead + 'timeseries_per_location_inclatlonadv_ddeg%d_sp%d_dd%d.nc'%(ddeg, sp, dd))
ts = pf[var][:]
#tssalin = pf['salin'][:]
tslon = pf['lon'][:]
tslat = pf['lat'][:]

Lons = pf['Lons'][:]
Lats = pf['Lats'][:]
vLons, vLats = np.meshgrid(Lons, Lats)
vLons = vLons.reshape(vLons.size)
vLats = vLats.reshape(vLats.size)
#% import the fixed averages (without advection)
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
    
dirReadfix = '/OFES/OFESres/TM/output/box-meandiff/'

pffix = np.load(dirReadfix +  'TM_box-meandiff'+var+'_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd))+'.npz')
trmfix = pffix['TMfix'][:]

zon_in_TM = region_index(lats.astype(int), lons.astype(int)).astype(np.int)
temp_nadv = trmfix[zon_in_TM]
ts = ts[zon_in_TM];

#Revise the order of the clusters from cold to warm
avgtemps = []
cl = range(n_clusters)
for i in range(n_clusters):
    avgtemps.append(np.mean(temp_nadv[czon==i]))
    newcl = [x for _,x in sorted(zip(avgtemps,cl))]
for i in range(n_clusters):
    print np.where(kmeans.labels_==newcl[i])
    print i
    czon[kmeans.labels_==newcl[i]] = i    

bd = {'temp ($^{\circ}$C)': temp_nadv[temp_nadv!=0],
      'cluster':czon[temp_nadv!=0] +1
}

#%% Final subplot
tail_perc = 5

fig = plt.figure(figsize=(20,10))
plt.subplots_adjust(hspace=0.33)

sns.set(style='darkgrid', font_scale=1.7)

plt.subplot2grid((2,2), (0, 0), colspan=1,rowspan=2)
plt.title('(a)', style='italic')
m = Basemap(projection= 'splaea',boundinglat=-12,lon_0=90,resolution='l')
m.drawcoastlines()
m.drawmapboundary(fill_color='w')
m.fillcontinents(color='grey')
m.drawparallels(np.arange(-70, 70, 20), labels=[True, False, False, False])
m.drawmeridians(np.arange(0, 361, 45), labels=[False, False, False, True]) 

xzon, yzon = m(lons, lats)

for i in range(n_clusters):
    idx = np.where(czon==i)
    plt.scatter(xzon[idx], yzon[idx], s=70, c=colors[i], zorder=3, label='cluster %d'%(i+1), marker=symbols[i])
plt.legend(title='Cluster', loc=2)

ax = plt.subplot2grid((2,2), (0, 1))

ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")

bspecies = [72]
plt.ylabel('mean LPDF($^{\circ}$C)')
cmap    = plt.get_cmap('nipy_spectral')

for i in range(len(bspecies)):
    idx = bspecies[i] - 19
    ra_species = sdata[:,idx]
    clusters = czon[ra_species==0]
    temp = temp_nadv[ra_species==0]   
    ra_species = ra_species[ra_species==0]
    ra_species = ra_species[temp>0]; clusters = clusters[temp>0]; temp = temp[temp>0]; 

    bd_species = {
    'Relative abundace (%)':  np.append(ra_species, np.zeros(1)),
    'cluster':np.append(clusters+1,np.array([1])),
    '$^{\circ}$C': np.append(temp,np.array([-50]))
    }

    ra_species = sdata[:,idx]
    clusters = czon[ra_species>0]
    temp = temp_nadv[ra_species>0]#[temp_nadv!=0]    
    ra_species = ra_species[ra_species>0]    
    bd_species2 = {
    'Relative abundace (%)':  ra_species,
    'cluster':clusters+1,
    '$^{\circ}$C': temp
    }    
    plt.title('(b) '+species[idx][:-3], style='italic')
    
    colo = {}#LogNorm  Normalize
    norm=matplotlib.colors.LogNorm(vmin = 0.1, vmax=100)
    for cval in bd_species2['Relative abundace (%)']:
        colo.update({cval : cmap(norm(cval))})
    
    sns.swarmplot(x='cluster',y='$^{\circ}$C', hue='Relative abundace (%)', size=9, data=bd_species, color='k')#, cmap='jet', palette=colors)#, whis=[0,100])
    sns.swarmplot(x='cluster',y='$^{\circ}$C', hue='Relative abundace (%)', size=9, data=bd_species2, palette=colo)#, whis=[0,100])
    plt.ylim(-3,30)
    plt.gca().legend_.remove()
    
    
    plt.xlabel('Cluster', x=0.1)
    
ax = plt.subplot2grid((2,2), (1, 1))
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
cluster = 1
species_no = 72 - 19
ra_species = sdata[:,species_no]

n_loc = np.sum(czon==cluster); idx_loc = np.where(czon==cluster);
left_tail = []
right_tail = []
mean_LPDF = []
ra = []
for i in range(n_loc):
    temp_series = ts[idx_loc[0][i]]
    temp_series = temp_series[temp_series<200]
    if(len(temp_series)>0):
        ra.append(ra_species[idx_loc[0][i]])
        left_tail.append(np.nanpercentile(temp_series,tail_perc))
        right_tail.append(np.nanpercentile(temp_series,100-tail_perc))
        mean_LPDF.append(temp_nadv[idx_loc[0][i]])

mean_LPDF= np.array(mean_LPDF); left_tail = np.array(left_tail); right_tail = np.array(right_tail); ra= np.array(ra)

print species[species_no], '   and cluster  ', cluster +1

id_ra = np.array(ra)>0
plt.scatter(mean_LPDF, np.array(left_tail) - np.array(mean_LPDF) , c='k', s=70)
plt.scatter(mean_LPDF[id_ra], np.array(left_tail)[id_ra] - np.array(mean_LPDF)[id_ra] , c=np.array(ra)[id_ra], cmap=cmap, s=70, vmin=0.1, vmax=100,  norm=matplotlib.colors.LogNorm())
plt.ylabel('AB1$\mathregular{^{cold}}$ ($^{\circ}$C)')
plt.title('(c) Cluster 2', style='italic')
plt.xlabel('mean LPDF ($^{\circ}$C)')

cbar = plt.colorbar(orientation='horizontal', cax = fig.add_axes([0.57, 0.01, 0.3, 0.03]),  label='Relative abundance(%)')

cbar.ax.plot([0.005,0.005],[0,1], 'k', linewidth=4.5)
plt.savefig('figure8.pdf', bbox_inches="tight")
plt.show()