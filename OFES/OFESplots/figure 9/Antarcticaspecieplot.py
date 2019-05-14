# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 13:30:56 2018

@author: nooteboom
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from netCDF4 import Dataset
import matplotlib.colors as colors
import Basemap

def running_mean(l, N):
    sum = 0
    result = list( 0 for x in l) 
    for i in range( 0, N ):
        sum = sum + l[i]
        result[i] = sum / (i+1) 
    for i in range( N, len(l) ):
        sum = sum - l[i-N] + l[i]
        result[i] = sum / N 
    return result

def running_mean_circular(l, N):
    leng = len(l)
    l = np.concatenate((l,l,l))
    sum = 0
    result = list( 0 for x in l) 
    for i in range( 0, N ):
        sum = sum + l[i]
        result[i] = sum / (i+1) 
    for i in range( N, len(l) ):
        sum = sum - l[i-N] + l[i]
        result[i] = sum / N 
    return result[leng:2*leng]

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

font = {'family' : 'Helvetica',
#        'weight' : 'bold',
        'size'   : 16}
ifont = {'family' : 'Helvetica',
        'style' : 'italic',
        'size'   : 18}
matplotlib.rc('font', **font) 

sp = 6
dd = 10
ddeg = 1

var = 'temp'

redcolor = 'firebrick'

savefile = False
delicatus = True
resolution = 'l'

#define functions
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

#%% Load Zonneveld data:
zonread = '/Users/nooteboom/Documents/PhD/parcels/secondpart/'

zondata = np.load(zonread + 'Zonneveldt2013_data.npy').item()

#%%Load the model PDFs after advection
dirRead = '/Users/nooteboom/Documents/PhD/parcels/OFES_global/Transtition_matrices/input/'

#pf = np.load(dirRead + 'timeseries_per_location_inclatlonadv_ddeg%d_sp%d_dd%d.npz'%(ddeg, sp, dd))
pf = Dataset(dirRead + 'timeseries_per_location_inclatlonadv_ddeg%d_sp%d_dd%d.nc'%(ddeg, sp, dd))
ts = pf[var][:]
tssalin = pf['salin'][:]
tslon = pf['lon'][:]
tslat = pf['lat'][:]

Lons = pf['Lons'][:]
Lats = pf['Lats'][:]
vLons, vLats = np.meshgrid(Lons, Lats)
vLons = vLons.reshape(vLons.size)
vLats = vLats.reshape(vLats.size)

#%% import the fixed averages (without advection)
dirReadfix = '/Users/nooteboom/Documents/PhD/parcels/OFES_global/Transtition_matrices/output/box-meandiff/'

pffix = np.load(dirReadfix +  'TM_box-meandiff'+var+'_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd))+'.npz')
trmfix = pffix['TMfix'][:]

pffixsalin = np.load(dirReadfix +  'TM_box-meandiff'+'salin'+'_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd))+'.npz')
trmfixsalin = pffixsalin['TMfix'][:]
#%% Locations of a sepcie
lb = 0.   #scatter only relative abundances above 'lb'

#Interesting cold species: 72, 40, 25

# 72: S. Antarctica

for specie_number in [72]:#[72]:#

    print 'specie: ', specie_number, '    ', zondata['variables'][specie_number]
    relative_abundances = zondata['data'][:,specie_number]
    sites_index = relative_abundances.astype(float) > 0
    relative_abundances = relative_abundances[sites_index].astype(float)
    
    zonlat = zondata['data'][:,1].astype(float)[sites_index]
    zonlon = zondata['data'][:,2].astype(float)[sites_index]; zonlon[np.where(zonlon<0)] += 360   
    
    zonlat = zonlat[relative_abundances>lb]
    zonlon = zonlon[relative_abundances>lb]
    relative_abundances = relative_abundances[relative_abundances>lb]
    
    zon_in_TM = region_index(zonlat, zonlon).astype(np.int) #indices where the zonneveld measurements are located in TM
    
    
ts = ts[zon_in_TM]
tssalin = tssalin[zon_in_TM]
tslon = tslon[zon_in_TM]
tslat = tslat[zon_in_TM]

temp_nadv = trmfix[zon_in_TM]
salin_nadv = trmfixsalin[zon_in_TM]

LI0 = np.zeros(ts.shape[0])
LI = np.zeros(ts.shape[0])
LI_salin = np.zeros(ts.shape[0])
mean = np.zeros(ts.shape[0])
ra = np.zeros(ts.shape[0])
glon = np.array([])
glat = np.array([])
notglon = np.array([])
notglat = np.array([])
for i in range(ts.shape[0]):
    if(len(ts[i])>0):
        ra[i] = relative_abundances[i] #relative abundance at measurement location i

        LI0[i] = np.percentile(np.array(ts[i]),ra[i])        
        LI[i] = np.percentile(np.array(ts[i]),ra[i]/2.) # nth percentile of PDF associated with measurement location i (Here n the relative abundance at measurement location i)     
        glon = np.append(glon, tslon[i][np.array(ts[i])< LI0[i]])
        glat = np.append(glat, tslat[i][np.array(ts[i])< LI0[i]])
        LI_salin[i] = np.nanmedian( tssalin[i][np.array(ts[i])< LI0[i]])
        notglon = np.append(notglon, tslon[i][np.array(ts[i])>= LI0[i]])
        notglat = np.append(notglat, tslat[i][np.array(ts[i])>= LI0[i]])
        

#%% Make the plot

fig, ax = plt.subplots(figsize=(12,8))

ACCloc = np.load('ACClocation.npz')
acclat = ACCloc['acclat'][:]
acclon = ACCloc['acclon'][:]
acclat = running_mean_circular(acclat,150)# To get a running_mean over circular data

plt.subplot2grid((16, 8), (0, 0), colspan=4, rowspan=12)

m = Basemap(projection='splaea',boundinglat=-30,lon_0=90,resolution=resolution)#north pole: 'nplaea', south pole: 'spaeqd'
m.drawcoastlines()
m.drawmapboundary(fill_color='w')#black')
m.fillcontinents(color='grey')#, lake_color='lightskyblue')
m.drawparallels(np.arange(-70, 70, 20), labels=[True, False, False, False])
m.drawmeridians(np.arange(0, 361, 45), labels=[False, False, False, True]) 
plt.title('(a) ' + zondata['variables'][specie_number][:-3],**ifont)

aclon, aclat = m(acclon, acclat)
xzon, yzon = m(zonlon, zonlat)
xadv, yadv = m(glon, glat)
xadv2, yadv2 = m(notglon, notglat)

cmap = plt.get_cmap('Blues')

plt.scatter(xadv, yadv, c=redcolor, label='backtracked cold tail', s=15, alpha=0.3)#, zorder=4

plt.plot(aclon,aclat, c='k', linewidth=2, label='ACC')#, zorder=2
plt.scatter(xzon, yzon, s=18, c=relative_abundances, cmap=truncate_colormap(cmap, 0.4, 1.0),   norm=colors.LogNorm(vmin=0.1, vmax=100))#vmin=-20, vmax=30, zorder = 0,
cb = plt.colorbar(orientation='horizontal', cax = fig.add_axes([0.145, 0.21, 0.34, 0.03]) )
cb.set_label('Measurement relative abundance(%)', labelpad=-60, rotation=0)#, y=1.05)

plt.scatter(xadv, yadv, c=redcolor, label='backtracked cold tail', s=15, alpha=0.3)#, zorder=4
plt.plot(aclon,aclat, c='k', linewidth=2, label='ACC')#, zorder=2
plt.legend(bbox_to_anchor=(0.8, -1.2, 0.1, .08))

plt.subplot2grid((16, 8), (0, 5), rowspan=6, colspan=3)
#scatter plot
plt.scatter(temp_nadv,relative_abundances, s=15, c='navy', label='mean LPDF', alpha=0.7)
plt.scatter(LI, ra, s=15, color=redcolor, zorder=10, label='tail APDF', alpha=0.7)
plt.ylabel('relative abundance')
plt.xlim(-3,23)
plt.title('(b) ' + zondata['variables'][specie_number],**ifont)
#%% Last subplot for equatorial specie
ts = pf[var][:]
tslon = pf['lon'][:]
tslat = pf['lat'][:]

if(delicatus):
    for specie_number in [74]:
    
        print 'specie: ', specie_number, '    ', zondata['variables'][specie_number]
        relative_abundances = zondata['data'][:,specie_number]
        sites_index = relative_abundances.astype(float) > 0
        relative_abundances = relative_abundances[sites_index].astype(float)
        
        zonlat = zondata['data'][:,1].astype(float)[sites_index]
        zonlon = zondata['data'][:,2].astype(float)[sites_index]; zonlon[np.where(zonlon<0)] += 360

        
        zon_in_TM = region_index(zonlat, zonlon).astype(np.int) #indices where the zonneveld measurements are located in TM
        
        
    ts = ts[zon_in_TM]
    tslon = tslon[zon_in_TM]
    tslat = tslat[zon_in_TM]
    temp_nadv = trmfix[zon_in_TM]
    
    LI = np.zeros(ts.shape[0])
    mean = np.zeros(ts.shape[0])
    ra = np.zeros(ts.shape[0])
    glon = np.array([])
    glat = np.array([])
    notglon = np.array([])
    notglat = np.array([])
    for i in range(ts.shape[0]):
        if(len(ts[i])>0):
            ra[i] = relative_abundances[i] #relative abundance at measurement location i
     
            LI[i] = np.percentile(np.array(ts[i]),100-ra[i]/2.) # 100-nth percentile of PDF associated with measurement location i (Here n the relative abundance at measurement location i)
            glon = np.append(glon, tslon[i][np.array(ts[i])<= LI[i]])
            glat = np.append(glat, tslat[i][np.array(ts[i])<= LI[i]])
            notglon = np.append(notglon, tslon[i][np.array(ts[i])>= LI[i]])
            notglat = np.append(notglat, tslat[i][np.array(ts[i])>= LI[i]])        
            
    lb = 0        

plt.subplot2grid((16, 8), (10, 5), rowspan=6, colspan=3)
if(delicatus):
    plt.scatter(temp_nadv,relative_abundances, s=15, c='navy', label='mean LPDF', alpha=0.7)
    plt.scatter(LI[ra>lb], ra[ra>lb], s=15, color=redcolor, zorder=10, label='tail APDF', alpha=0.7)
    plt.xlim(3,33)
    plt.xlabel('temperature ($^{\circ}C$)')
    plt.title('(c) ' + zondata['variables'][specie_number],**ifont)    
else:
    #scatterplot of salinity:
    plt.scatter(salin_nadv,relative_abundances, s=15, c='navy', label='mean LPDF', alpha=0.7)
    plt.scatter(LI_salin[LI_salin<100], ra[LI_salin<100], s=15, color=redcolor, zorder=10, label='tail APDF', alpha=0.7)
    plt.xlabel('salinity (PSU)')
    plt.title('(c) ' ,**ifont)

plt.ylabel('relative abundance')

plt.legend(bbox_to_anchor=(0.8, 1.49, 0.1, .08))

if(savefile):
    plt.savefig('tailpdf_forspecie.pdf', bbox_inches="tight")
plt.show()
#%%
trywarmspecie = True
#Interesting for equator: [87, 83, 81, 74, 73, 49, 33]
        # especially: 73, 74, 33
if(trywarmspecie):
    ts = pf[var][:]
    tslon = pf['lon'][:]
    tslat = pf['lat'][:]
    
    lonconstrain = True
    if(lonconstrain):
        lowerlon = 360-65
        upperlon = 360-0
        lowerlat = 30

    for specie_number in [74]:
    
        print 'specie: ', specie_number, '    ', zondata['variables'][specie_number]
        relative_abundances = zondata['data'][:,specie_number]
        sites_index = relative_abundances.astype(float) > 0
        relative_abundances = relative_abundances[sites_index].astype(float)
        
        zonlat = zondata['data'][:,1].astype(float)[sites_index]
        zonlon = zondata['data'][:,2].astype(float)[sites_index]; zonlon[np.where(zonlon<0)] += 360
        zontemp = zondata['data'][:,8].astype(float)[sites_index]
        
        if(lonconstrain):
            index = np.where(np.logical_and(np.logical_and(zonlon<=upperlon, zonlon>lowerlon),zonlat>=lowerlat))
            zonlat = zonlat[index]
            zonlon = zonlon[index]
            zontemp = zontemp[index]            
            relative_abundances = relative_abundances[index]
            if(upperlon>360):
                zonlon[zonlon<100] += 360
        
        zon_in_TM = region_index(zonlat, zonlon).astype(np.int) #indices where the zonneveld measurements are located in TM
        
        
    ts = ts[zon_in_TM]
    tslon = tslon[zon_in_TM]
    tslat = tslat[zon_in_TM]
    
    temp_nadv = trmfix[zon_in_TM]
    
    LI = np.zeros(ts.shape[0])
    mean = np.zeros(ts.shape[0])
    ra = np.zeros(ts.shape[0])
    glon = np.array([])
    glat = np.array([])
    notglon = np.array([])
    notglat = np.array([])       

    for i in range(ts.shape[0]):
        if(len(ts[i])>0):
            ra[i] = relative_abundances[i] #relative abundance at measurement location i
    #        print ra[i]      
            LI[i] = np.percentile(np.array(ts[i]),100-ra[i]/2.) # 100-nth percentile of PDF associated with measurement location i (Here n the relative abundance at measurement location i)
            glon = np.append(glon, tslon[i][np.array(ts[i])>= LI[i]])
            glat = np.append(glat, tslat[i][np.array(ts[i])>= LI[i]])
            notglon = np.append(notglon, tslon[i][np.array(ts[i])<= LI[i]])
            notglat = np.append(notglat, tslat[i][np.array(ts[i])<= LI[i]])   

            if( lonconstrain):
                if(upperlon>360):
                    glon[glon<100] += 360
                    notglon[notglon<100] += 360        
    
    #% Make the plot
    
    plt.subplots(figsize=(40,12))
    
    plt.subplot(121)
    
    if(lonconstrain):
        m = Basemap(projection='cyl', llcrnrlat=0, urcrnrlat=60, llcrnrlon=lowerlon-5, urcrnrlon=upperlon+10, resolution=resolution)  
    else:
        m = Basemap(projection='cyl', llcrnrlat=-75, urcrnrlat=75, llcrnrlon=0, urcrnrlon=360, resolution=resolution) 
    m.drawcoastlines()
    m.drawmapboundary(fill_color='w')#black')
    m.fillcontinents(color='grey')#, lake_color='lightskyblue')
    m.drawparallels(np.arange(-90, 90, 60), labels=[True, False, False, False])
    m.drawmeridians(np.arange(0, 361, 90), labels=[False, False, False, True]) 
    
    xzon, yzon = m(zonlon, zonlat)
    xadv, yadv = m(glon, glat)
    xadv2, yadv2 = m(notglon, notglat)
    
    plt.scatter(xadv, yadv, s=15, alpha=0.3, c = 'b')
    plt.scatter(xzon, yzon, c=temp_nadv, s=15, vmin=13, vmax=17, cmap = 'spectral')
    plt.colorbar(orientation = 'horizontal')
    
    plt.subplot(122)
    #scatter plot
    lb = 0    #scatter only relative abundances above 'lb'
    
    plt.scatter(temp_nadv,relative_abundances, s=15, c='blue', label='mean LPDF')
    plt.scatter(LI[ra>lb], ra[ra>lb], s=15, color=redcolor, zorder=10, label='tail APDF')
    
    plt.xlabel('temperature ($^{\circ}C$)')
    plt.ylabel('relative abundance (%)')
    plt.xlim(0.2,40)
    plt.legend()
    
    plt.show()