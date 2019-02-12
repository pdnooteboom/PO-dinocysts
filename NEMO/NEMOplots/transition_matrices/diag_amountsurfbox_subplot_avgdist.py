# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 15:22:53 2018

Plot the diagonal of the transition matrix and the amount of boxes any bottom box is mapped to.
Use 6 m/s here.
 

@author: nooteboom
"""

import numpy as np
import matplotlib.pylab as plt
import matplotlib
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmocean
import matplotlib.colors as colors

plotcmap = 'pink'#'plasma'#cmocean.cm.cmap_d[]
plotcmap2 = 'bone'#'viridis'#cmocean.cm.cmap_d['haline'] # 
plotcmap3 = 'copper'#'winter'#'viridis'#'viridis'#cmocean.cm.cmap_d['haline'] # 
cmap3 = plt.get_cmap(plotcmap3)
plotresolution = 'l'

# True if surface instead of amount of surface grid boxes
surfacebool = True

ddeg = 2 # resolution of the binning
sp = 6
dd = 10
res = 1

tmdir = '/Users/nooteboom/Documents/PhD/parcels/NEMO/atsf/Transition_matrices/'

data = np.load(tmdir + 'output/box-box/TMglobal_bin'+str(ddeg)+'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res) + '.npz')
TM = data['TM'][:]
Lons = data['Lons'][:]
Lats = data['Lats'][:]
#%%Load avg drift distance
datadist = np.load(tmdir + 'output/box-avgdist/TM_box-avgdist__ddeg%d_sp%d_dd%d.npz'%(ddeg, sp, dd))
TMd = datadist['TM'][:]
Lonsd = datadist['Lons'][:]
Latsd = datadist['Lats'][:]
vLons, vLats = np.meshgrid(Lonsd, Latsd)
vLons = vLons.flatten(); vLats = vLats.flatten();

#Lonsd[Lonsd>=180] -= 360; 

#TMd = TMd[np.logical_and(vLons>-0.5,vLats<=86)]
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
#%% Use only part of colormap
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
        
#%%

font = {'family' : 'Helvetica',
'weight' : 'normal',
'size'   : 20}

matplotlib.rc('font', **font) 

class oceanvector(object):

    
    def __init__(self,vec,vec2, vec3, Lons=[],Lats=[], Lons2=[],Lats2=[],val=None):
        
        if not np.any(Lons):
            self.Lons = np.linspace(280.,325., 325-280)
            self.Lats = np.linspace(25., 47., 47-25) 
        else:
            self.Lons = Lons
            self.Lats = Lats
            self.Lons2 = Lons2
            self.Lats2 = Lats2
            
        if vec.ndim==1:
            v1d = vec
            v2d = vec.reshape((len(self.Lats),len(self.Lons)))
        else:
            v1d = vec.ravel()
            v2d = vec
            
        if vec.ndim==1:
            v1d2 = vec2
            v2d2 = vec2.reshape((len(self.Lats),len(self.Lons)))
        else:
            v1d2 = vec2.ravel()
            v2d2 = vec2          

        if vec.ndim==1:
            v1d3 = vec3
            v2d3 = vec3.reshape((len(self.Lats2),len(self.Lons2)))
#            v2d3 = np.concatenate((v2d3[:,91:],v2d3[:,:91]), axis=1)
        else:
            v1d3 = vec3.ravel()
            v2d3 = vec3  
        
        self.V1d = v1d
        self.V2d = v2d  
        self.V1d2 = v1d2
        self.V2d2 = v2d2        
        self.V1d3 = v1d3
        self.V2d3 = v2d3 
        
    def plot_me(self, bounds = False, vbounds = None, vbounds3 = None, land= True, cbartitle='', cbartitle2='', cbartitle3='' ,
                cmap='inferno', cmap2='viridis', cmap3='viridis', colbar = True, colbar2=True, orien='vertical', title = None, 
                pl='cmesh', save = False, outname=''):
                    
        parallelplots = [-75,-25,25,70]
        meridianplots = [-120,0,120]
        
        if(save):
            respl = 'l'
        else:
            respl = 'l'
            
        matplotlib.rc('font', **font)        
        fig = plt.figure(figsize=(15,20))
#First subplot
        ax = fig.add_subplot(311)
        plt.title('(a)')
        if(bounds):
            m = Basemap(projection='cyl', llcrnrlat=bounds[1][0], urcrnrlat=bounds[1][1], llcrnrlon=bounds[0][0], urcrnrlon=bounds[0][1], resolution=respl)
        else:
            m = Basemap(projection='cyl', llcrnrlat=np.min(self.Lats), urcrnrlat=np.max(self.Lats), llcrnrlon=np.min(self.Lons), urcrnrlon=np.max(self.Lons), resolution=plotresolution)
        m.drawcoastlines()
        m.drawparallels(parallelplots,labels=[True,False,True,False])#np.arange(np.min(self.Lats),np.max(self.Lats),50.)
        m.drawmeridians(meridianplots,labels=[False,False,False,False])
        if land:
            m.fillcontinents(color='silver')
            m.drawmapboundary(fill_color='k')
        lon_bins_2d, lat_bins_2d = np.meshgrid(self.Lons, self.Lats)
    
        xs, ys = m(lon_bins_2d, lat_bins_2d)
            
        #Large halo for plotting for lon>360
        xs = np.concatenate((xs,xs+358,xs+2*358), axis=1)
        ys = np.concatenate((ys,ys,ys), axis=1) 
        zs = self.V2d
        zs = np.concatenate((zs,zs,zs), axis=1)            
        
        if pl=='cmesh':                   
            plt.pcolormesh(xs, ys, zs*100, cmap=cmap, vmin=0., vmax=100.) #,cmap=cm.viridis)
        else:
            plt.contour(xs, ys, zs,cmap='set1') #,cmap=cm2.coolwarm)
        
        if colbar:
            if orien=='vertical':
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cbar = plt.colorbar(orientation=orien, cax=cax)
            else:
                cbaxes = fig.add_axes([0.125, 0.27, 0.352, 0.03])
                cbar = plt.colorbar(orientation=orien, cax=cbaxes)
            cbar.set_label(cbartitle,size=20) 

#Second subplot            
        ax2 = fig.add_subplot(312) 
        plt.title('(b)')         
        if(bounds):
            m = Basemap(projection='cyl', llcrnrlat=bounds[1][0], urcrnrlat=bounds[1][1], llcrnrlon=bounds[0][0], urcrnrlon=bounds[0][1], resolution=respl)
        else:
            m = Basemap(projection='cyl', llcrnrlat=np.min(self.Lats), urcrnrlat=np.max(self.Lats), llcrnrlon=np.min(self.Lons), urcrnrlon=np.max(self.Lons), resolution=plotresolution)
        m.drawcoastlines()

        m.drawparallels(parallelplots,labels=[True,False,False,False])#np.arange(np.min(self.Lats),np.max(self.Lats)+1,50.)
        m.drawmeridians(meridianplots,labels=[False,False,False,False])
        if land:
            m.fillcontinents(color='silver')
            m.drawmapboundary(fill_color='k')
        lon_bins_2d, lat_bins_2d = np.meshgrid(self.Lons, self.Lats)
     
        zs = self.V2d2
        zs = np.concatenate((zs,zs,zs), axis=1)            
        
        if pl=='cmesh':                   
            plt.pcolormesh(xs, ys, zs, cmap=cmap2, vmin=vbounds[0], vmax=vbounds[1])#norm=colors.LogNorm(vmin=vbounds[0], vmax=vbounds[1]))#
        else:
            plt.contour(xs, ys, zs,cmap='set1') #,cmap=cm2.coolwarm)
               
        if colbar2:
            if orien=='vertical':
                divider = make_axes_locatable(ax2)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cbar = plt.colorbar(orientation=orien, cax=cax)
            else:
                cbaxes = fig.add_axes([0.548, 0.27, 0.352, 0.03])
                cbar = plt.colorbar(orientation=orien, cax=cbaxes)
            cbar.set_label(cbartitle2,size=20)

#Third subplot            
        ax3 = fig.add_subplot(313) 
        plt.title('(c)')         
        if(bounds):
            m = Basemap(projection='cyl', llcrnrlat=bounds[1][0], urcrnrlat=bounds[1][1], llcrnrlon=bounds[0][0], urcrnrlon=bounds[0][1], resolution=respl)
        else:
            m = Basemap(projection='cyl', llcrnrlat=np.min(self.Lats2), urcrnrlat=np.max(self.Lats2), llcrnrlon=np.min(self.Lons2), urcrnrlon=np.max(self.Lons2), resolution=plotresolution)
        m.drawcoastlines()

        m.drawparallels(parallelplots,labels=[True,False,False,False])#np.arange(np.min(self.Lats),np.max(self.Lats)+1,50.)
        m.drawmeridians(meridianplots,labels=[False,False,False,True])
        if land:
            m.fillcontinents(color='silver')
            m.drawmapboundary(fill_color='k')
        lon_bins_2d, lat_bins_2d = np.meshgrid(self.Lons2, self.Lats2)
        
        xs, ys = m(lon_bins_2d, lat_bins_2d)  
        xs = np.concatenate((xs,xs+358), axis=1)
        ys = np.concatenate((ys,ys), axis=1)         
        zs = self.V2d3
        zs = np.concatenate((zs,zs), axis=1)            
        
        if pl=='cmesh':                   
            plt.pcolormesh(xs, ys, zs, cmap=cmap3, vmin=vbounds3[0], vmax=vbounds3[1])#, norm=colors.LogNorm(vmin=vbounds3[0], vmax=vbounds3[1]))##
        else:
            plt.contour(xs, ys, zs,cmap='set1') #,cmap=cm2.coolwarm)
               
        if colbar2:
            if orien=='vertical':
                divider = make_axes_locatable(ax3)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cbar = plt.colorbar(orientation=orien, cax=cax)
            else:
                cbaxes = fig.add_axes([0.548, 0.27, 0.352, 0.03])
                cbar = plt.colorbar(orientation=orien, cax=cbaxes)
            cbar.set_label(cbartitle3,size=20)
                   
#General           
         
        if title is None:
            print 'no title'
        else:
            plt.title(title,size=18)            
            
        if save:
            plt.savefig(outname,bbox_inches='tight')
            plt.close()
        else:
            plt.show()     
            

#%%

bounds = [[180,540],[-75,70]]        
threshold = 0.
connected = (TM>threshold).astype(int)
print connected.shape
if(surfacebool):
    for i in range(connected.shape[0]):
        connected[i] = connected[i] * surface
    cbartitle2 = '$10^6\ km^2$'
    vbounds = [0, 8] 
else:
    cbartitle2 = '(#)'
    vbounds = [0, 150]  
    


N = oceanvector(np.diagonal(TM, offset=0),np.sum(connected, axis=0)/1000000., TMd/1000.,Lons=Lons,Lats=Lats,Lons2=Lonsd,Lats2=Latsd)
np.save('TM_contour', N)
cbartitle3 = '$10^3\ km$'
vbounds3 = [0.1,2.5]

N.plot_me(bounds = bounds, vbounds = vbounds, vbounds3 = vbounds3, land= True, cbartitle='(%)', cbartitle2=cbartitle2, cbartitle3 = cbartitle3, 
                cmap=plotcmap, cmap2 = plotcmap2, cmap3 = truncate_colormap(cmap3, 0.1, 1.), colbar = True, colbar2=True, orien='vertical', title = None, 
                pl='cmesh', save =True, outname = '/Users/nooteboom/Documents/PhD/firstpaper/articleplots/nemoplots/' + 'TM_NEMO_diag_amountboxes_avgdist.eps')        
        