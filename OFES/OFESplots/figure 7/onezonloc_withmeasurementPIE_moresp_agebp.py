# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 17:05:38 2018

@author: nooteboom
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import numpy.ma as ma
import seaborn as sns
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties

font = {'family' : 'Helvetica',
        'weight' : 'bold',
        'size'   : 18}

matplotlib.rc('font', **font) 

def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

resolution = 'l'
savefile = True
if(savefile):
    resolution = 'l'
    
redcolor = 'firebrick'
greencolor = 'mediumseagreen'#'seagreen'
yellowcolor = 'goldenrod'#'gold'

sns.set(color_codes=True)
sns.set_context("poster")
sns.set(font_scale=1.4)

#%%
exdir = '/OFES/OFESres/wholetraj/'

sp = 6
sp2=11
sp3=50

#aus:
plon = 154.91
plat =  -49.71
dd = 10

outfile = '2loc_dd10_sp'+ str(int(sp)) + '_lon'+ str(plon) + '_lat'+ str(plat) + ".nc"
outfile_sp2 = 'loc_dd10_sp'+ str(int(sp2)) + '_lon'+ str(plon) + '_lat'+ str(plat) + ".nc"
outfile_sp3 = 'loc_dd10_sp'+ str(int(sp3)) + '_lon'+ str(plon) + '_lat'+ str(plat) + ".nc"

pfile2 = Dataset(exdir  + outfile, 'r')  
pfile2_sp2 = Dataset(exdir  + outfile_sp2, 'r')
pfile2_sp3 = Dataset(exdir + outfile_sp3, 'r')
    
time = np.unique(pfile2['time'])[::-1]

#%% Open Zonneveld data:
zonread = ''

zondata = np.load(zonread + 'Zonneveldt2013_data.npy').item()
zonlat = zondata['data'][:,1].astype(float)
zonlon = zondata['data'][:,2].astype(float)
reg = np.where(zonlat==plat)#np.where(np.logical_and(zonlat==plat,zonlon+360==plon))
rel_abundances =  zondata['data'][reg[0],19:-2][0].astype(float)
species = np.array(zondata['variables'][19:-2])[np.where(rel_abundances>0)]
rel_abundances = np.array(rel_abundances)[np.where(rel_abundances>0)[0]]

print 'total amount of species: ', len(species)
#%% preprocessing before plot
day = np.arange(0,3*len(pfile2['time'][0,:]),3)

tempind = np.array([], dtype = int)

slon = np.array([])
slat = np.array([])
stemp = np.array([])
sage = np.array([])

tempind_sp2 = np.array([], dtype = int)

slon_sp2 = np.array([])
slat_sp2 = np.array([])
stemp_sp2 = np.array([])
sage_sp2 = np.array([])

tempind_sp3 = np.array([], dtype = int)

slon_sp3 = np.array([])
slat_sp3 = np.array([])
stemp_sp3 = np.array([])
sage_sp3 = np.array([])

for t in range(len(day)):

    ind = np.where(np.logical_and(pfile2['z'][:,t]>10.,ma.getmask(pfile2['z'][:,-1])))
    lon = pfile2['lon'][:,t][ind]
    lat = pfile2['lat'][:,t][ind]
    z = pfile2['z'][:,t][ind]
       
    tempind = np.where(np.invert(ma.getmaskarray(pfile2['temp'][:,t])))[0].astype(int)
    
    if(tempind.shape[0]>0):
        slon = np.append(slon,pfile2['lon'][:,t][(tempind)])
        slat = np.append(slat,pfile2['lat'][:,t][(tempind)])
        
        stemp = np.append(stemp,pfile2['temp'][:,t][(tempind)])    
        sage = np.append(sage,pfile2['age'][:,t][(tempind)])    

    ends_ind = ma.getmask(pfile2['z'][:,-1])    
    
    if(pfile2_sp2['temp'][:].shape[1]>t):
        tempind_sp2 = np.where(np.invert(ma.getmaskarray(pfile2_sp2['temp'][:,t])))[0].astype(int)

        if(tempind_sp2.shape[0]>0):
            slon_sp2 = np.append(slon_sp2,pfile2_sp2['lon'][:,t][(tempind_sp2)])
            slat_sp2 = np.append(slat_sp2,pfile2_sp2['lat'][:,t][(tempind_sp2)])
            
            stemp_sp2 = np.append(stemp_sp2,pfile2_sp2['temp'][:,t][(tempind_sp2)])
            sage_sp2 = np.append(sage_sp2,pfile2_sp2['age'][:,t][(tempind_sp2)])   

        ends_ind_sp2 = ma.getmask(pfile2_sp2['z'][:,-1])

    if(pfile2_sp3['temp'][:].shape[1]>t):
        tempind_sp3 = np.where(np.invert(ma.getmaskarray(pfile2_sp3['temp'][:,t])))[0].astype(int)

        if(tempind_sp3.shape[0]>0):
            slon_sp3 = np.append(slon_sp3,pfile2_sp3['lon'][:,t][(tempind_sp3)])
            slat_sp3 = np.append(slat_sp3,pfile2_sp3['lat'][:,t][(tempind_sp3)])
            
            stemp_sp3 = np.append(stemp_sp3,pfile2_sp3['temp'][:,t][(tempind_sp3)])
            sage_sp3 = np.append(sage_sp3,pfile2_sp3['age'][:,t][(tempind_sp3)])   
        ends_ind_sp3 = ma.getmask(pfile2_sp3['z'][:,-1])
#%% plot general
font = {'family' : 'Helvetica',
'weight' : 'bold',
'size'   : 15}

matplotlib.rc('font', **font)

fig = plt.figure(figsize=(17,10))
# set up subplot grid
gridspec.GridSpec(18, 18)
#%%First subplot in lat-lon plane

plt.subplot2grid((18,18), (0,0), colspan=8, rowspan=6)
m = Basemap(projection='cyl', llcrnrlat=-70, urcrnrlat=-35, llcrnrlon=90, urcrnrlon=170, resolution=resolution) 
m.drawcoastlines()
m.drawmapboundary(fill_color='w')
m.fillcontinents(color='grey')
m.drawparallels(np.arange(-70, 70, 10), labels=[True, False, False, False])
m.drawmeridians(np.arange(0, 350, 20), labels=[False, False, False, True]) 

xs, ys = slon, slat
xs0, ys0 = [plon], [plat] 
xs_sp2, ys_sp2 = slon_sp2, slat_sp2
xs_sp3, ys_sp3 = slon_sp3, slat_sp3

print '# particles: ', len(xs),len(xs_sp2), len(xs_sp3)

plt.scatter(xs, ys, c=redcolor, s=40, alpha=0.28)
plt.scatter(xs_sp2, ys_sp2, c=greencolor, s=40, alpha=0.28)
plt.scatter(xs_sp3, ys_sp3, c=yellowcolor, s=40, alpha=0.28)
plt.scatter(xs0, ys0, c='navy', s=180, marker = 'P')

plt.title('(a)')
#%%
font = {'family' : 'Helvetica',
#'weight' : 'bold',
'size'   : 10}

matplotlib.rc('font', **font)
sage = sage/3600./24./365.
sage_sp2 = sage_sp2/3600./24/365.
sage_sp3 = sage_sp3/3600./24./365.

print '% above 1 year travel time sp6:   ', np.sum(sage>1)/np.float(len(sage))
print '% above 1 year travel time sp11:   ',np.sum(sage_sp2>1)/np.float(len(sage_sp2))
data = { 'sinking speed':np.concatenate((np.array([sp]*len(sage[sage<4])),np.array([sp2]*len(sage_sp2)),np.array([sp3]*len(sage_sp3)))),
        'Time (years)':np.concatenate((sage[sage<4],sage_sp2, sage_sp3))
        
        }

plt.subplot2grid((18,18), (7,0), colspan=8, rowspan=11)
plt.title('(c)')
ax = sns.violinplot(y='sinking speed', x='Time (years)', 
                    palette={sp:redcolor, sp2:greencolor, sp3:yellowcolor}, data=data, orient = 'h', 
                    inner='quartile', scale = 'count',
                    width = 0.95)
plt.xlabel('Travel time (years)')
plt.ylabel('(m day$^{-1}$)')
#%% relative abundances of species
font = {'family' : 'Helvetica',
'weight' : 'bold',
'size'   : 12}

matplotlib.rc('font', **font)
fontP = FontProperties()
fontP.set_size('medium')

#% Pie chart
plt.subplot2grid((18,18), (9,10), colspan=5, rowspan=8)

colormap = plt.cm.gist_rainbow   

size = 0.8
vals = rel_abundances
outer_colors = [colormap(i) for i in np.linspace(0, 0.9,len(species))]

plt.pie(vals, radius=1, colors=outer_colors,startangle=80,
       wedgeprops=dict(width=size, edgecolor='k'), autopct='%1.1f%%',pctdistance=1.2, center=(0,3))#, center=(-5,0)

plt.legend(np.array([spe[:-4] for spe in species]),loc='center',ncol=2, bbox_to_anchor=(0.5, -0.2), prop=fontP)#{'size': 6})
plt.title('(d)\n')
#%% Third subplot: pdfs of temperature (could also add salinity etc.)
font = {'family' : 'Helvetica',
'weight' : 'bold',
'size'   : 15}

matplotlib.rc('font', **font)
sns.set_style("darkgrid")
plt.subplot2grid((18,18), (0,9), colspan=8, rowspan=6)

pfile = Dataset(exdir + 'surfaceloc_lon'+ str(plon) + '_lat'+ str(plat) + '.nc')
    
stemp_nadv = ma.getdata(pfile['temp'][:])
stemp_nadv = stemp_nadv[0,1:]

print 'mean of pdf advected (6 m/day): ',np.nanmean(stemp)
print '5-th and 95-th percentile around median:   ', np.nanpercentile(stemp, 5), np.nanpercentile(stemp,95)
print '5-th and 95-th percentile around median (no advection):   ', np.nanpercentile(stemp_nadv, 5), np.nanpercentile(stemp_nadv,95)
print 'mean of pdf advected (11 m/day): ',np.nanmean(stemp_sp2)
print 'mean of pdf advected (50 m/day): ',np.nanmean(stemp_sp3)
print 'mean of pdf fixed: ',np.nanmean(stemp_nadv)
print 'median of pdf fixed: ',np.nanmedian(stemp_nadv)

plt.axvline(x=np.nanmean(stemp), dashes=[6,2], color=redcolor, linewidth= 3)
plt.axvline(x=np.nanmean(stemp_sp2), dashes=[6,2], color=greencolor, linewidth= 3)
plt.axvline(x=np.nanmean(stemp_sp3), dashes=[6,2], color=yellowcolor, linewidth= 3)
plt.axvline(x=np.nanmean(stemp_nadv), dashes=[6,2], color='navy', linewidth= 3)

sns.kdeplot(stemp, bw=1.5, linewidth = 5, label='6 m day$^{-1}$', color=redcolor);
sns.kdeplot(stemp_sp2, bw=1.5, linewidth = 5, label='11 m day$^{-1}$', color=greencolor);
sns.kdeplot(stemp_sp3, bw=1.5, linewidth = 5, label='50 m day$^{-1}$', color=yellowcolor);
sns.kdeplot(stemp_nadv, bw=1.5, linewidth = 5, label='Local', color='navy');

print 'standard deiations of advected temperature and fixed temperature: ', np.std(stemp), '    ',np.std(stemp_nadv)

plt.xlabel('Sea Surface Temperature ($^{\circ}$ C)', fontsize='19', horizontalalignment='center')#, fontweight='bold'
plt.xticks( fontweight='bold', fontsize='17', horizontalalignment='right')
plt.yticks( fontweight='bold', fontsize='17', horizontalalignment='right')
plt.legend(loc='b', fontsize=17)
plt.title('(b)')
plt.xlim(-8, 17)
if(savefile):
    plt.savefig('station3627.pdf', bbox_inches="tight")
plt.show()