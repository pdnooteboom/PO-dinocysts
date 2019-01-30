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
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import seaborn as sns


plotcmap = 'pink'
plotcmap2 = 'bone'
plotcmap3 = 'copper'
cmap3 = plt.get_cmap(plotcmap3)
plotresolution = 'l'

# True if surface instead of amount of surface grid boxes
surfacebool = True
dd = 100
sp = 6#'s2'
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

class SeabornFig2Grid():

    def __init__(self, seaborngrid, fig,  subplot_spec):
        self.fig = fig
        self.sg = seaborngrid
        self.subplot = subplot_spec
        if isinstance(self.sg, sns.axisgrid.FacetGrid) or \
            isinstance(self.sg, sns.axisgrid.PairGrid):
            self._movegrid()
        elif isinstance(self.sg, sns.axisgrid.JointGrid):
            self._movejointgrid()
        self._finalize()

    def _movegrid(self):
        """ Move PairGrid or Facetgrid """
        self._resize()
        n = self.sg.axes.shape[0]
        m = self.sg.axes.shape[1]
        self.subgrid = gridspec.GridSpecFromSubplotSpec(n,m, subplot_spec=self.subplot)
        for i in range(n):
            for j in range(m):
                self._moveaxes(self.sg.axes[i,j], self.subgrid[i,j])

    def _movejointgrid(self):
        """ Move Jointgrid """
        h= self.sg.ax_joint.get_position().height
        h2= self.sg.ax_marg_x.get_position().height
        r = int(np.round(h/h2))
        self._resize()
        self.subgrid = gridspec.GridSpecFromSubplotSpec(r+1,r+1, subplot_spec=self.subplot)

        self._moveaxes(self.sg.ax_joint, self.subgrid[1:, :-1])
        self._moveaxes(self.sg.ax_marg_x, self.subgrid[0, :-1])
        self._moveaxes(self.sg.ax_marg_y, self.subgrid[1:, -1])

    def _moveaxes(self, ax, gs):
        #https://stackoverflow.com/a/46906599/4124317
        ax.remove()
        ax.figure=self.fig
        self.fig.axes.append(ax)
        self.fig.add_axes(ax)
        ax._subplotspec = gs
        ax.set_position(gs.get_position(self.fig))
        ax.set_subplotspec(gs)

    def _finalize(self):
        plt.close(self.sg.fig)
        self.fig.canvas.mpl_connect("resize_event", self._resize)
        self.fig.canvas.draw()

    def _resize(self, evt=None):
        self.sg.fig.set_size_inches(self.fig.get_size_inches())


class oceanvector(object):

    
    def __init__(self,vec,vec2, vec3, depth, Lons=[],Lats=[],val=None):
        
        if not np.any(Lons):
            self.Lons = np.linspace(280.,325., 325-280)
            self.Lats = np.linspace(25., 47., 47-25) 
        else:
            self.Lons = Lons
            self.Lats = Lats
        
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
            v2d3 = vec3.reshape((len(self.Lats),len(self.Lons)))
        else:
            v1d3 = vec3.ravel()
            v2d3 = vec3  
        
        self.V1d = v1d
        self.V2d = v2d  
        self.V1d2 = v1d2
        self.V2d2 = v2d2        
        self.V1d3 = v1d3
        self.V2d3 = v2d3 
        self.depth = depth
        
    def plot_me(self, bounds = False, vbounds = None, vbounds3 = None, ub=90, land= True, cbartitle='', cbartitle2='', cbartitle3='' ,
                cmap='inferno', cmap2='viridis', cmap3='viridis', colbar = True, colbar2=True, orien='vertical', title = None, 
                pl='cmesh', save = False, outname=''):
          
        lon_bins_2d, lat_bins_2d = np.meshgrid(self.Lons, self.Lats)
        scat_lat1 = lat_bins_2d.flatten()<ub        
    
        
        if(save):
            respl = 'l'
        else:
            respl = 'l'
            
        matplotlib.rc('font', **font)        
        fig = plt.figure(figsize=(18,18))
# set up subplot grid
        gs = gridspec.GridSpec(3, 3,
                               width_ratios=[140,1, 30],
                               height_ratios=[1, 1, 1]
                               )
#First subplot        
        ax = plt.subplot(gs[0])
        plt.title('(a)')
        if(bounds):
            m = Basemap(projection='cyl', llcrnrlat=bounds[1][0], urcrnrlat=bounds[1][1], llcrnrlon=bounds[0][0], urcrnrlon=bounds[0][1], resolution=respl)
        else:
            m = Basemap(projection='cyl', llcrnrlat=np.min(self.Lats), urcrnrlat=np.max(self.Lats), llcrnrlon=np.min(self.Lons), urcrnrlon=np.max(self.Lons), resolution=plotresolution)
        m.drawcoastlines()
        m.drawparallels([np.min(self.Lats),np.min(self.Lats)+50,np.min(self.Lats)+100,70],labels=[True,False,True,False])#np.arange(np.min(self.Lats),np.max(self.Lats),50.)
        m.drawmeridians(np.arange(np.min(self.Lons),np.max(self.Lons),120.),labels=[False,False,False,False])
        if land:
            m.fillcontinents(color='silver')
            m.drawmapboundary(fill_color='k')
        lon_bins_2d, lat_bins_2d = np.meshgrid(self.Lons, self.Lats)
    
        xs, ys = m(lon_bins_2d, lat_bins_2d)
            
        #Large halo for plotting for lon>360
        xs = np.concatenate((xs,xs+358), axis=1)
        ys = np.concatenate((ys,ys), axis=1)        
        zs1 = self.V2d
        zs = np.concatenate((zs1,zs1), axis=1)            
        
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
            cbar.ax.set_title(cbartitle,size=20) 

#Second subplot            
        ax2 = plt.subplot(gs[3])
        plt.title('(b)')         
        if(bounds):
            m = Basemap(projection='cyl', llcrnrlat=bounds[1][0], urcrnrlat=bounds[1][1], llcrnrlon=bounds[0][0], urcrnrlon=bounds[0][1], resolution=respl)
        else:
            m = Basemap(projection='cyl', llcrnrlat=np.min(self.Lats), urcrnrlat=np.max(self.Lats), llcrnrlon=np.min(self.Lons), urcrnrlon=np.max(self.Lons), resolution=plotresolution)
        m.drawcoastlines()

        m.drawparallels([np.min(self.Lats),np.min(self.Lats)+50,np.min(self.Lats)+100,70],labels=[True,False,False,False])#np.arange(np.min(self.Lats),np.max(self.Lats)+1,50.)
        m.drawmeridians(np.arange(np.min(self.Lons),np.max(self.Lons),120.),labels=[False,False,False,False])
        if land:
            m.fillcontinents(color='silver')
            m.drawmapboundary(fill_color='k')
        lon_bins_2d, lat_bins_2d = np.meshgrid(self.Lons, self.Lats)
     
        zs2 = self.V2d2

        zs = np.concatenate((zs2,zs2), axis=1)            
        
        if pl=='cmesh':                   
            plt.pcolormesh(xs, ys, zs, cmap=cmap2, vmin=vbounds[0], vmax=vbounds[1])#norm=colors.LogNorm(vmin=vbounds[0], vmax=vbounds[1]))#
        else:
            plt.contour(xs, ys, zs,cmap='set1') #,cmap=cm2.coolwarm)
               
        if colbar2:
            if orien=='vertical':
                divider = make_axes_locatable(ax2)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cbar = plt.colorbar(orientation=orien, cax=cax)#, location = 1.0)
            else:
                cbaxes = fig.add_axes([0.548, 0.27, 0.352, 0.03])
                cbar = plt.colorbar(orientation=orien, cax=cbaxes)
            cbar.ax.set_title(cbartitle2,size=20)

#Third subplot            
        ax3 = plt.subplot(gs[6])        
        plt.title('(c)')         
        if(bounds):
            m = Basemap(projection='cyl', llcrnrlat=bounds[1][0], urcrnrlat=bounds[1][1], llcrnrlon=bounds[0][0], urcrnrlon=bounds[0][1], resolution=respl)
        else:
            m = Basemap(projection='cyl', llcrnrlat=np.min(self.Lats), urcrnrlat=np.max(self.Lats), llcrnrlon=np.min(self.Lons), urcrnrlon=np.max(self.Lons), resolution=plotresolution)
        m.drawcoastlines()

        m.drawparallels([np.min(self.Lats),np.min(self.Lats)+50,np.min(self.Lats)+100,70],labels=[True,False,False,False])#np.arange(np.min(self.Lats),np.max(self.Lats)+1,50.)
        m.drawmeridians(np.arange(np.min(self.Lons),np.max(self.Lons),120.),labels=[False,False,False,True])
        if land:
            m.fillcontinents(color='silver')
            m.drawmapboundary(fill_color='k')
        lon_bins_2d, lat_bins_2d = np.meshgrid(self.Lons, self.Lats)
     
        zs3 = self.V2d3
        zs = np.concatenate((zs3,zs3), axis=1)            
        
        if pl=='cmesh':                   
            plt.pcolormesh(xs, ys, zs, cmap=cmap3, vmin=vbounds3[0], vmax=vbounds3[1])#, norm=colors.LogNorm(vmin=vbounds3[0], vmax=vbounds3[1]))##
        else:
            plt.contour(xs, ys, zs,cmap='set1') 
               
        if colbar2:
            if orien=='vertical':
                divider = make_axes_locatable(ax3)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cbar = plt.colorbar(orientation=orien, cax=cax)
            else:
                cbaxes = fig.add_axes([0.548, 0.27, 0.352, 0.03])
                cbar = plt.colorbar(orientation=orien, cax=cbaxes)
            cbar.ax.set_title(cbartitle3,size=20)

        depths = self.depth

        # scatter plots
            
        ssize = 3
        ec='k'
        aa = 0.2  
        yl = (-0.2,6.2)
        kin = 'scatter'#'hex'#
               

        if(kin=='scatter'):        
            g = sns.jointplot(x=zs1.flatten()[scat_lat1]*100, y= depths[scat_lat1]/1000., kind=kin, color=ec, s=ssize, alpha=aa, ylim = yl, xlim = (0,102))#, ax=ax)
        else:        
            g = sns.jointplot(x=zs1.flatten()[scat_lat1]*100, y= depths[scat_lat1]/1000., kind=kin, color=ec, ylim = yl, xlim = (0,102))#, ax=ax)

        g.fig.suptitle('(d)')
        plt.gca().invert_yaxis()
        g.ax_joint.set_xlabel(cbartitle)#, fontweight='bold')
        g.ax_joint.set_ylabel('depth (km)')#, fontweight='bold')
        SeabornFig2Grid(g, fig, gs[2])

        if(kin=='scatter'):        
            g = sns.jointplot(x=zs2.flatten()[scat_lat1], y= depths[scat_lat1]/1000., kind=kin, color=ec, s=ssize, alpha=aa, ylim = yl, xlim = (0,9))
        else:        
            g = sns.jointplot(x=zs2.flatten()[scat_lat1], y= depths[scat_lat1]/1000., kind=kin, color=ec, ylim = yl, xlim = (0,9))
        plt.gca().invert_yaxis()
        g.ax_joint.set_xlabel(cbartitle2)#, fontweight='bold')
        g.ax_joint.set_ylabel('depth (km)')#, fontweight='bold')
        SeabornFig2Grid(g, fig, gs[5])

        if(kin=='scatter'):
            g = sns.jointplot(x=zs3.flatten()[scat_lat1], y= depths[scat_lat1]/1000., kind=kin, color=ec, s=ssize, alpha=aa, ylim = yl, xlim = (0,2.8))#, ax=ax)
        else:
            g = sns.jointplot(x=zs3.flatten()[scat_lat1], y= depths[scat_lat1]/1000., kind=kin, color=ec, ylim = yl, xlim = (0,2.8))#, ax=ax)


            
        plt.gca().invert_yaxis()
        g.ax_joint.set_xlabel(cbartitle3)#, fontweight='bold')
        g.ax_joint.set_ylabel('depth (km)')#, fontweight='bold')

        SeabornFig2Grid(g, fig, gs[8])

  
        gs.tight_layout(fig)

      
#General           
         
        if title is None:
            print 'no title'
        else:
            plt.title(title,size=18)            
            
        if save:
            print 'fig dpi:  ',fig.dpi
#            fig.savefig(outname[:-3]+'pdf', dpi=fig.dpi, bbox_inches='tight')#, pad_inches=15)#,bbox_inches='tight')
            fig.savefig(outname[:-3]+'png', dpi=fig.dpi, bbox_inches='tight')#, pad_inches=15)#,bbox_inches='tight')            
            plt.show()
#            plt.close()
        else:
            plt.show()     
            

#%%


    
ub = 71 #upper latitude bound for the scatter plot

bounds = [[180,540],[-75,70]]   
cbartitle='(%)'
cbartitle2 = '$10^6\ km^2$'    
cbartitle3 = '$10^3\ km$'

pd = np.load('process/toplot_TM_diag_amountboxes_avgdist_scatter_sp'+str(sp)+'_dd'+str(dd) +'.npz')
Lons = pd['Lons']
Lats = pd['Lats']  
N = oceanvector(pd['TM'],pd['connected']/1000000., pd['TMd']/1000., pd['depth'],Lons=Lons,Lats=Lats)
outname =  'TM_diag_amountboxes_avgdist_jointplot_sp'+str(sp)+'_dd%d2.png'%(int(dd))
vbounds = [0, 8] 
vbounds3 = [0.,2.5]

N.plot_me(bounds = bounds, vbounds = vbounds, vbounds3 = vbounds3, ub=ub, land= True, cbartitle=cbartitle, cbartitle2=cbartitle2, cbartitle3 = cbartitle3, 
                cmap=plotcmap, cmap2 = plotcmap2, cmap3 = truncate_colormap(cmap3, 0.1, 1.), colbar = True, colbar2=True, orien='vertical', title = None, 
                pl='cmesh', save =True, outname = outname)        
        