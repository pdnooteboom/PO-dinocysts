# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 15:22:53 2018

Plot the diagonal of the transition matrix and the amount of boxes any bottom box is mapped to.
Use 6 m/s here.
 

@author: nooteboom
"""

import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib
import seaborn as sns

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

    
    def __init__(self,vec,vec2, vec3, depth, depthd, Lons=[],Lats=[], Lons2=[],Lats2=[],val=None):
        
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
        self.depth = depth.reshape((len(self.Lats),len(self.Lons)))
        self.depthd = depthd.reshape((len(self.Lats2),len(self.Lons2)))
        
    def plot_me(self, bounds = False, vbounds = None, vbounds3 = None, ub=90, land= True, cbartitle='', cbartitle2='', cbartitle3='' ,
                cmap='inferno', cmap2='viridis', cmap3='viridis', colbar = True, colbar2=True, orien='vertical', title = None, 
                pl='cmesh', save = False, outname=''):

        ssize = 3
        ec='darkblue'
        aa = 0.1
                    
        parallelplots = [-75,-25,25,70]
        meridianplots = [-120,0,120]
        
        if(save):
            respl = 'l'
        else:
            respl = 'l'
            
        matplotlib.rc('font', **font)    
        
        fig = plt.figure(figsize=(18,18))
# set up subplot grid
#        gs = gridspec.GridSpec(18, 24)
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
            m = Basemap(projection='cyl', llcrnrlat=np.min(self.Lats), urcrnrlat=np.max(self.Lats), llcrnrlon=np.min(self.Lons), urcrnrlon=np.max(self.Lons), resolution='l')
        m.drawcoastlines()
        m.drawparallels(parallelplots,labels=[True,False,True,False])#np.arange(np.min(self.Lats),np.max(self.Lats),50.)
        m.drawmeridians(meridianplots,labels=[False,False,False,False])
        if land:
            m.fillcontinents(color='silver')
            m.drawmapboundary(fill_color='k')
        lon_bins_2d, lat_bins_2d = np.meshgrid(self.Lons, self.Lats)

        scat_lat1 = lat_bins_2d.flatten()<ub
    
        xs, ys = m(lon_bins_2d, lat_bins_2d)
            
        #Large halo for plotting for lon>360
        xs = np.concatenate((xs,xs+358,xs+2*358), axis=1)
        ys = np.concatenate((ys,ys,ys), axis=1) 
        zs1 = self.V2d
        zs = np.concatenate((zs1,zs1,zs1), axis=1) 
        depth1 = self.depth
#        depth = np.concatenate((depth1,depth1,depth1), axis=1)   
        
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
        ax2 = plt.subplot(gs[3])
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
     
        zs2 = self.V2d2
        zs = np.concatenate((zs2,zs2,zs2), axis=1)            
        
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
        ax3 = plt.subplot(gs[6])  
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
        
        scat_lat2 = lat_bins_2d.flatten()<ub
        
        xs, ys = m(lon_bins_2d, lat_bins_2d)  
        xs = np.concatenate((xs,xs+358), axis=1)
        ys = np.concatenate((ys,ys), axis=1)         
        zs3 = self.V2d3
        zs = np.concatenate((zs3,zs3), axis=1)
        depthd3 = self.depthd
#        depthd = np.concatenate((depthd3,depthd3), axis=1)             
        
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

        # scatter plots
            
        ssize = 3
        ec='royalblue'#'viridis'#'k'
        aa = 0.2   
        yl = (-0.2,6.2)
        kin = 'scatter'#'kde'#'hex'#scatter
        #pal = sns.cubehelix_palette(light=1, as_cmap=True)#sns.dark_palette("palegreen", as_cmap=True)
        pal = sns.set_palette(reversed(sns.color_palette("Blues_d", 20)), 20)
               
#        ax = plt.subplot2grid((18,24), (1,19), colspan=6, rowspan=4)
#        self.plot_joint(title = '(d)\n', cbartitle=cbartitle, zs = zs1.flatten()*100, depths= depths, subplotspec=gs[1])

#        g = sns.JointGrid(x=zs1.flatten()*100, y= depths/1000.)#, ax=ax)
#        g = g.plot_joint(plt.scatter, color="white", edgecolor="k")
#        g = g.plot_marginals(sns.distplot, kde=False, color=".5")
        if(kin=='scatter'):
            g = sns.jointplot(x=zs1.flatten()[scat_lat1]*100, y= depth1.flatten()[scat_lat1]/1000., kind=kin, color=ec, s=ssize, alpha=aa, ylim = yl, xlim = (0,102))
        else:
            sns.set(style="white")
            g = sns.jointplot(x=zs1.flatten()[scat_lat1]*100, y= depth1.flatten()[scat_lat1]/1000., kind=kin, color=ec, ylim = yl, xlim = (0,102), cmap = pal)
        #plt.title('(d)')
        g.fig.suptitle('(d)')
        plt.gca().invert_yaxis()
        g.ax_joint.set_xlabel(cbartitle)#, fontweight='bold')
        g.ax_joint.set_ylabel('depth (km)')#, fontweight='bold')
        SeabornFig2Grid(g, fig, gs[2])

        if(kin=='scatter'):
            g = sns.jointplot(x=zs2.flatten()[scat_lat1], y= depth1.flatten()[scat_lat1]/1000., kind=kin, color=ec, s=ssize, alpha=aa, ylim = yl, xlim = (0,9))
        else:
            g = sns.jointplot(x=zs2.flatten()[scat_lat1], y= depth1.flatten()[scat_lat1]/1000., kind=kin, color=ec, ylim = yl, xlim = (0,9), cmap = pal)
        plt.gca().invert_yaxis()
        g.ax_joint.set_xlabel(cbartitle2)#, fontweight='bold')
        g.ax_joint.set_ylabel('depth (km)')#, fontweight='bold')
        SeabornFig2Grid(g, fig, gs[5])

        if(kin=='scatter'):
            g = sns.jointplot(x=zs3.flatten()[scat_lat2], y= depthd3.flatten()[scat_lat2]/1000., kind=kin, color=ec, s=ssize, alpha=aa, ylim = yl, xlim = (0,2.8))
        else:
            g = sns.jointplot(x=zs3.flatten()[scat_lat2], y= depthd3.flatten()[scat_lat2]/1000., kind=kin, color=ec, ylim = yl, xlim = (0,2.8), cmap = pal)
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
            plt.savefig(outname)#,bbox_inches='tight')
            plt.show()
        else:
            plt.show()     
            

#%%
plotcmap = 'pink'#'plasma'#cmocean.cm.cmap_d[]
plotcmap2 = 'bone'#'viridis'#cmocean.cm.cmap_d['haline'] # 
plotcmap3 = 'copper'#'winter'#'viridis'#'viridis'#cmocean.cm.cmap_d['haline'] # 
cmap3 = plt.get_cmap(plotcmap3)
plotresolution = 'l'

surfacebool = True

bounds = [[180,540],[-75,70]]        
if(surfacebool):
    cbartitle2 = '$10^6\ km^2$'
    vbounds = [0, 8] 
else:
    cbartitle2 = '(#)'
    vbounds = [0, 150]  

pd = np.load('process/toplot_TM_NEMO_diag_amountboxes_avgdist_scatter.npz')
Lons = pd['Lons']
Lats = pd['Lats']
Lonsd = pd['Lons2']
Latsd = pd['Lats2']
ub = 71 #upper latitude bound for the scatter plot

N = oceanvector(pd['TM'],pd['connected']/1000000., pd['TMd']/1000., pd['depth'], pd['depthd'],Lons=Lons,Lats=Lats,Lons2=Lonsd,Lats2=Latsd)
cbartitle3 = '$10^3\ km$'
vbounds3 = [0.1,2.5]



N.plot_me(bounds = bounds, vbounds = vbounds, vbounds3 = vbounds3, ub=ub, land= True, cbartitle='(%)', cbartitle2=cbartitle2, cbartitle3 = cbartitle3, 
                cmap=plotcmap, cmap2 = plotcmap2, cmap3 = truncate_colormap(cmap3, 0.1, 1.), colbar = True, colbar2=True, orien='vertical', title = None, 
                pl='cmesh', save =False, outname = '/Users/nooteboom/Documents/PhD/firstpaper/articleplots/nemoplots/' + 'TM_NEMO_diag_amountboxes_avgdist_joint_hex.png')        
        