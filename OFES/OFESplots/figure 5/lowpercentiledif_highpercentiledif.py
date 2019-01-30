# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 10:27:36 2018

@author: nooteboom
"""

import numpy as np
import matplotlib.pylab as plt
import matplotlib
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
import cmocean
import scipy.ndimage as ndimage

plotres = 'l'

plotcontour = True
contourrange = [-2,2]
contourrange2 = [-0.5,0.5]
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
zoomam = 6
sigma = 1

ddeg = 1 # resolution of the binning at the bottom
sp = 6#'s1'#6
dd = 10
n = 5

pltse = 'annual'

var = 'temp'
var2 = 'salin'


plotcmap1 = cmocean.cm.cmap_d['balance']#['thermal'] # cmocean.cm.cmap_d['curl'] # 'Spectral_r'#'seismic'#
barmin1=-16; barmax1 = 16;barrange1 = range(barmin1,barmax1+1,4);

plotcmap2 = cmocean.cm.cmap_d['delta']#['haline'] # cmocean.cm.cmap_d['curl'] # 'Spectral_r'#'seismic'#
barmin2=-4; barmax2 = 4;barrange2 = range(barmin2,barmax2+1,2);

dirRead = '/Users/nooteboom/Documents/GitHub/PO-dinocysts/' + 'OFES/OFESres/TM/'


def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
    
def find_down(array,value):
    if(value>array[0]):
        array = [n-value for n in array]
        idx = np.array([n for n in array if n<0]).argmax()
    else:
        idx = 0
    return idx  
    
    
font = {'family' : 'Helvetica',
    'size'   : 22}

matplotlib.rc('font', **font) 

#%% Define figure
fig = plt.figure(figsize=(17,9))
fig.subplots_adjust(bottom=0.2)

#%% nth percentile var1
if(type(sp)==str):
    pf = np.load(dirRead +  'box-percentiledif/TM_box-mean%dpercentiledif_'%(n)+var+'_ddeg%d_sp'%(ddeg)+sp+'_dd%d'%(int(dd))+'.npz')
else:
    pf = np.load(dirRead +  'box-percentiledif/TM_box-mean%dpercentiledif_'%(n)+var+'_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd))+'.npz')

Lons = pf['Lons'][:]
Lats = pf['Lats'][:]

was = pf['TM'][:]
print was.shape, Lons.shape, Lats.shape

longitudes, latitudes = np.meshgrid(Lons, Lats)

was = was.reshape(longitudes.shape)
minplotlat = -75
maxplotlat = 70
minplotlon = 180
maxplotlon = 360+180
if(maxplotlon>360):
    minlati = find_down(latitudes[:,0], minplotlat)
    maxlati = find_down(latitudes[:,0], maxplotlat)
    minloni = find_down(longitudes[0], minplotlon)
    maxloni = find_down(longitudes[0]+360, maxplotlon)
    longitudes = np.concatenate((longitudes[:,maxloni:], longitudes[:,:maxloni]+360),axis=1)
    was = np.concatenate((was[:,maxloni:], was[:,:maxloni]),axis=1)
    minloni = find_down(longitudes[0], minplotlon)
    maxloni = find_down(longitudes[0], maxplotlon)
    
else:
    minlati = find_down(latitudes[:,0], minplotlat)
    maxlati = find_down(latitudes[:,0], maxplotlat)
    minloni = find_down(longitudes[0], minplotlon)
    maxloni = find_down(longitudes[0], maxplotlon)

ax = plt.subplot(221)
plt.title('(a)', fontsize=20)
m = Basemap(projection='cyl', llcrnrlat=latitudes[minlati,0], urcrnrlat=latitudes[maxlati,0], llcrnrlon=longitudes[0,minloni], urcrnrlon=longitudes[0,maxloni], resolution=plotres)
m.drawcoastlines(zorder=4)
m.fillcontinents(color='silver',zorder=3)
m.drawparallels(np.arange(-75, 76, 50), labels=[True, False, False, False],zorder=3)
m.drawmeridians(np.arange(-120, 180, 120), labels=[False, False, False, False],zorder=3) 

xs, ys = m(longitudes[minlati:maxlati+1, minloni:maxloni+1], latitudes[minlati:maxlati+1, minloni:maxloni+1])
zs = was[minlati:maxlati+1, minloni:maxloni+1]#[minloni:maxloni+1, minlati:maxlati+1]#
zs = np.ma.masked_where(np.isnan(zs),zs)
zs = zs[:-1,1:];xs = xs[:-1,1:]; ys = ys[:-1,1:];

x = xs[~zs.mask]
y = ys[~zs.mask]

nzs = griddata((x,y), zs[~zs.mask].ravel(), (xs,ys), method='linear')

m.imshow(nzs, extent=[np.min(xs), np.max(xs), np.min(ys), np.max(ys)], interpolation='spline36', cmap=plotcmap1,vmin=barmin1, vmax=barmax1,zorder=1)
if(plotcontour):
    CRt0 = plt.contour(xs,ys,ndimage.gaussian_filter(nzs, sigma=sigma, order=0),contourrange, colors='k',zorder=2)
#%%  nth percentile var2

if(type(sp)==str):
    pf = np.load(dirRead +  'box-percentiledif/TM_box-mean%dpercentiledif_'%(n)+var2+'_ddeg%d_sp'%(ddeg)+sp+'_dd%d'%(int(dd))+'.npz')   
else:
    pf = np.load(dirRead +  'box-percentiledif/TM_box-mean%dpercentiledif_'%(n)+var2+'_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd))+'.npz')

Lons = pf['Lons'][:]
Lats = pf['Lats'][:]

was = pf['TM'][:]
print was.shape, Lons.shape, Lats.shape

longitudes, latitudes = np.meshgrid(Lons, Lats)

was = was.reshape(longitudes.shape)
if(maxplotlon>360):
    minlati = find_down(latitudes[:,0], minplotlat)
    maxlati = find_down(latitudes[:,0], maxplotlat)
    minloni = find_down(longitudes[0], minplotlon)
    maxloni = find_down(longitudes[0]+360, maxplotlon)
    longitudes = np.concatenate((longitudes[:,maxloni:], longitudes[:,:maxloni]+360),axis=1)
    was = np.concatenate((was[:,maxloni:], was[:,:maxloni]),axis=1)
    minloni = find_down(longitudes[0], minplotlon)
    maxloni = find_down(longitudes[0], maxplotlon)    
else:
    minlati = find_down(latitudes[:,0], minplotlat)
    maxlati = find_down(latitudes[:,0], maxplotlat)
    minloni = find_down(longitudes[0], minplotlon)
    maxloni = find_down(longitudes[0], maxplotlon)

ax = plt.subplot(222)
plt.title('(c)', fontsize=20)
m = Basemap(projection='cyl', llcrnrlat=latitudes[minlati,0], urcrnrlat=latitudes[maxlati,0], llcrnrlon=longitudes[0,minloni], urcrnrlon=longitudes[0,maxloni], resolution=plotres)#
m.drawcoastlines(zorder=3)
m.fillcontinents(color='silver',zorder=2)
m.drawparallels(np.arange(-75, 76, 50), labels=[False, False, False, False])
m.drawmeridians(np.arange(-120, 180, 120), labels=[False, False, False, False]) 

xs, ys = m(longitudes[minlati:maxlati+1, minloni:maxloni+1], latitudes[minlati:maxlati+1, minloni:maxloni+1])
zs = was[minlati:maxlati+1, minloni:maxloni+1]#[minloni:maxloni+1, minlati:maxlati+1]#
zs = np.ma.masked_where(np.isnan(zs),zs)
zs = zs[:-1,1:];xs = xs[:-1,1:]; ys = ys[:-1,1:];

x = xs[~zs.mask]
y = ys[~zs.mask]

nzs = griddata((x,y), zs[~zs.mask].ravel(), (xs,ys), method='linear')

m.imshow(nzs, extent=[np.min(xs), np.max(xs), np.min(ys), np.max(ys)], interpolation='spline36', cmap=plotcmap2,vmin=barmin2, vmax=barmax2,zorder=1)
if(plotcontour):
    CRs0 = plt.contour(xs,ys,ndimage.gaussian_filter(nzs, sigma=sigma, order=0),contourrange2, colors='k',zorder=2)

#%% (100-n)th percentile var1

if(type(sp)==str):
    pf = np.load(dirRead +  'box-percentiledif/TM_box-mean%dpercentiledif_'%(100-n)+var+'_ddeg%d_sp'%(ddeg)+sp+'_dd%d'%(int(dd))+'.npz')
else:
    pf = np.load(dirRead +  'box-percentiledif/TM_box-mean%dpercentiledif_'%(100-n)+var+'_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd))+'.npz')

Lons = pf['Lons'][:]
Lats = pf['Lats'][:]

was = pf['TM'][:]
print was.shape, Lons.shape, Lats.shape

longitudes, latitudes = np.meshgrid(Lons, Lats)

was = was.reshape(longitudes.shape)
if(maxplotlon>360):
    minlati = find_down(latitudes[:,0], minplotlat)
    maxlati = find_down(latitudes[:,0], maxplotlat)
    minloni = find_down(longitudes[0], minplotlon)
    maxloni = find_down(longitudes[0]+360, maxplotlon)
    longitudes = np.concatenate((longitudes[:,maxloni:], longitudes[:,:maxloni]+360),axis=1)
    was = np.concatenate((was[:,maxloni:], was[:,:maxloni]),axis=1)
    minloni = find_down(longitudes[0], minplotlon)
    maxloni = find_down(longitudes[0], maxplotlon)
    
else:
    minlati = find_down(latitudes[:,0], minplotlat)
    maxlati = find_down(latitudes[:,0], maxplotlat)
    minloni = find_down(longitudes[0], minplotlon)
    maxloni = find_down(longitudes[0], maxplotlon)

ax = plt.subplot(223)
plt.title('(b)', fontsize=20)
m = Basemap(projection='cyl', llcrnrlat=latitudes[minlati,0], urcrnrlat=latitudes[maxlati,0], llcrnrlon=longitudes[0,minloni], urcrnrlon=longitudes[0,maxloni], resolution=plotres)
m.drawcoastlines(zorder=3)
m.fillcontinents(color='silver',zorder=2)
m.drawparallels(np.arange(-75, 76, 50), labels=[True, False, False, False])
m.drawmeridians(np.arange(-120, 180, 120), labels=[False, False, False, True]) 

xs, ys = m(longitudes[minlati:maxlati+1, minloni:maxloni+1], latitudes[minlati:maxlati+1, minloni:maxloni+1])
zs = was[minlati:maxlati+1, minloni:maxloni+1]
zs = np.ma.masked_where(np.isnan(zs),zs)
zs = zs[:-1,1:];xs = xs[:-1,1:]; ys = ys[:-1,1:];

x = xs[~zs.mask]
y = ys[~zs.mask]


nzs = griddata((x,y), zs[~zs.mask].ravel(), (xs,ys), method='linear')


im = m.imshow(nzs, extent=[np.min(xs), np.max(xs), np.min(ys), np.max(ys)], interpolation='spline36', cmap=plotcmap1,vmin=barmin1, vmax=barmax1,zorder=1)
if(plotcontour):
    plt.contour(xs,ys,ndimage.gaussian_filter(nzs, sigma=sigma, order=0),contourrange, colors='k',zorder=2)

cbar_ax = fig.add_axes([0.14, 0.1, 0.33, 0.05])
cbar = fig.colorbar(im, cax=cbar_ax, orientation = 'horizontal',ticks=barrange1)
cbar.ax.set_xlabel('$^{\circ}$C', fontdict={'size':25})
cbar.ax.tick_params(labelsize=22)
cbar.ax.plot([16/32.+2/32.]*2,[0,1], 'k')
cbar.ax.plot([16/32.-2/32.]*2,[0,1], 'k')
#%% (100-n)-th percentile of var2

if(type(sp)==str):
    pf = np.load(dirRead +  'box-percentiledif/TM_box-mean%dpercentiledif_'%(100-n)+var2+'_ddeg%d_sp'%(ddeg)+sp+'_dd%d'%(int(dd))+'.npz')
else:
    pf = np.load(dirRead +  'box-percentiledif/TM_box-mean%dpercentiledif_'%(100-n)+var2+'_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd))+'.npz')

Lons = pf['Lons'][:]
Lats = pf['Lats'][:]

was = pf['TM'][:]
print was.shape, Lons.shape, Lats.shape

longitudes, latitudes = np.meshgrid(Lons, Lats)

was = was.reshape(longitudes.shape)
if(maxplotlon>360):
    minlati = find_down(latitudes[:,0], minplotlat)
    maxlati = find_down(latitudes[:,0], maxplotlat)
    minloni = find_down(longitudes[0], minplotlon)
    maxloni = find_down(longitudes[0]+360, maxplotlon)
    longitudes = np.concatenate((longitudes[:,maxloni:], longitudes[:,:maxloni]+360),axis=1)
    was = np.concatenate((was[:,maxloni:], was[:,:maxloni]),axis=1)
    minloni = find_down(longitudes[0], minplotlon)
    maxloni = find_down(longitudes[0], maxplotlon)    
else:
    minlati = find_down(latitudes[:,0], minplotlat)
    maxlati = find_down(latitudes[:,0], maxplotlat)
    minloni = find_down(longitudes[0], minplotlon)
    maxloni = find_down(longitudes[0], maxplotlon)

ax = plt.subplot(224)
plt.title('(d)', fontsize=20)
m = Basemap(projection='cyl', llcrnrlat=latitudes[minlati,0], urcrnrlat=latitudes[maxlati,0], llcrnrlon=longitudes[0,minloni], urcrnrlon=longitudes[0,maxloni], resolution=plotres)
m.drawcoastlines(zorder=3)
m.fillcontinents(color='silver',zorder=2)
m.drawparallels(np.arange(-75, 76, 50), labels=[False, False, False, False])
m.drawmeridians(np.arange(-120, 180, 120), labels=[False, False, False, True]) 

xs, ys = m(longitudes[minlati:maxlati+1, minloni:maxloni+1], latitudes[minlati:maxlati+1, minloni:maxloni+1])
zs = was[minlati:maxlati+1, minloni:maxloni+1]
zs = np.ma.masked_where(np.isnan(zs),zs)
zs = zs[:-1,1:];xs = xs[:-1,1:]; ys = ys[:-1,1:];

x = xs[~zs.mask]
y = ys[~zs.mask]

nzs = griddata((x,y), zs[~zs.mask].ravel(), (xs,ys), method='linear')

im = m.imshow(nzs, extent=[np.min(xs), np.max(xs), np.min(ys), np.max(ys)], interpolation='spline36', cmap=plotcmap2,vmin=barmin2, vmax=barmax2,zorder=1)
if(plotcontour):
    plt.contour(xs,ys,ndimage.gaussian_filter(nzs, sigma=sigma, order=0),contourrange2, colors='k',zorder=2)

cbar_ax = fig.add_axes([0.56, 0.1, 0.33, 0.05])
cbar = fig.colorbar(im, cax=cbar_ax, orientation = 'horizontal',ticks=barrange2)
cbar.ax.set_xlabel('psu', fontdict={'size':25})
cbar.ax.tick_params(labelsize=22)
cbar.ax.plot([4/8.+0.5/8.]*2,[0,1], 'k')
cbar.ax.plot([4/8.-0.5/8.]*2,[0,1], 'k')

#%%
if(type(sp)==str):
    plt.savefig('dif_mean_low-highpercentile_sp'+sp+'_dd%d.pdf'%(int(dd)), bbox_inches="tight")
else:
    plt.savefig('dif_mean_low-highpercentile_sp'+str(sp)+'_dd%d.pdf'%(int(dd)), bbox_inches="tight")

plt.show()