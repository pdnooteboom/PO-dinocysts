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
import seaborn as sns
import pandas as pd
from matplotlib.font_manager import FontProperties

#%%
font = {'family' : 'Helvetica',
#'weight' : 'bold',
'size'   : 60}

matplotlib.rc('font', **font)

ddeg = 2 # resolution of the binning
spl = [3, 6, 's1','s2', 11, 25, 50, 100, 200, 500]
dd = 10
res = 1

#Set for Pdensity plots:
scale = 'exponential'
k_depth = 'proportion'
outlier_prop = 0.5
whis = [5, 95]

meanpointprops = dict(marker='D', markeredgecolor='black',
                      markerfacecolor='firebrick', markersize=7)

#%% boxplot of different sinking velocities

toplot = np.load('process/toplot_TM-boxplot.npz')


d = {'$10^6\ km^2$':toplot['areas'],
 'sinking speed (m day$^{-1}$)':toplot['sinkingspeed'],
' ':toplot['locations'] }   

df = pd.DataFrame(data=d)
sns.set(style="whitegrid", font='Helvetica', font_scale=1.2)

fig = plt.figure(figsize=(15,10))

ax2 = plt.subplot2grid((2,19),(0,10), colspan=9)


plt.title('(b)', fontsize=25)
order = [str(i) for i in spl];
for i in range(len(order)): 
    if(order[i]=='s1'): 
        order[i] = 'SC1'
    elif(order[i]=='s2'):
        order[i] = 'SC2'
order = np.array(order)
ax = sns.boxplot(x='sinking speed (m day$^{-1}$)', y='$10^6\ km^2$', hue=' ', data=df, showmeans=True, meanprops=meanpointprops, order=order, whis=whis,  showfliers=True)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles[:2], labels=labels[:2])

#% Boxplot of the diagonal of the transition matrix


d =     {'%':toplot['diagonal'],
         'sinking speed (m day$^{-1}$)': toplot['sinkingspeed_diag'],
        ' ': toplot['locations_diag']}
        
np.save('TM_boxplot_firstsubplot',d)

df = pd.DataFrame(data=d)

sns.set(style="whitegrid", font='Helvetica', font_scale=1.2)
ax2 = plt.subplot2grid((2,19),(0,0), colspan=9)
ax = sns.boxplot(x='sinking speed (m day$^{-1}$)', y='%', hue=' ', data=df, showmeans=True, meanprops=meanpointprops, order=order, whis=whis,  showfliers=False)
ax.legend_.remove()
plt.title('(a)', fontsize=25)



d =     {'km':toplot['distance'],
         'sinking speed (m day$^{-1}$)': toplot['sinkingspeed_dist'],
        ' ': toplot['locations_dist']}

df = pd.DataFrame(data=d)        
ax2 = plt.subplot2grid((2,19),(1,5), colspan=9)
ax = sns.boxplot(x='sinking speed (m day$^{-1}$)', y='km', hue=' ', data=df, showmeans=True, meanprops=meanpointprops, order=order, whis=whis,  showfliers=False)
ax.legend_.remove()
plt.title('(c)', fontsize=25)       


plt.title('(d)', fontsize=25)    

ax3 = plt.subplot2grid((2,38),(0,0), colspan=12)#9)
plt.title('(a)', fontsize=25)  
plt.xlabel('sinking speed (m day$^{-1}$)       ')

res = 10000

z = np.linspace(10,6000,res)
c3 = np.full(res, 3)
c6 = np.full(res, 6)
c11 = np.full(res, 11)
c25 = np.full(res, 25)
c50 = np.full(res, 50)
c100 = np.full(res, 100)
c200 = np.full(res, 200)
c500 = np.full(res, 500)

def Sink1(z):
    if(z<100):
        return 6
    elif(z<2000):
        return 6 + (45-6)*(z-100)/np.float(2000-100)
    else:
        return 45
        
def Sink2(z):
    if(z<100):
        return 6
    elif(z<2000):
        return 6 + (45-6)*(z-100)/np.float(2000-100)
    elif(z<3500):
        return 45 + (65-45)*(z-2000)/np.float(3500-2000)
    else:
        return 65
        
cSC1 = np.array([Sink1(i) for i in z])
cSC2 = np.array([Sink2(i) for i in z])

linewidth = 3

z = z / 1000.

plt.plot(cSC1,z, linewidth = linewidth, label = 'SC1')
plt.plot(cSC2,z,'--', linewidth = linewidth, label = 'SC2')

plt.plot(c3,z, linewidth = linewidth, label = '3m day$^{-1}$')
plt.plot(c6,z, linewidth = linewidth, label = '6m day$^{-1}$')
plt.plot(c11,z, linewidth = linewidth, label = '11m day$^{-1}$')
plt.plot(c25,z, linewidth = linewidth, label = '25m day$^{-1}$')
plt.plot(c50,z, linewidth = linewidth, label = '50m day$^{-1}$')
plt.plot(c100,z, linewidth = linewidth, label = '100m day$^{-1}$')
plt.plot(c200,z, linewidth = linewidth, label = '200m day$^{-1}$')
plt.plot(c500,z, linewidth = linewidth, label = '500m day$^{-1}$')


fontP = FontProperties()
fontP.set_size('small')
plt.legend(bbox_to_anchor=(1.01, 1.0))#,['SC1', 'SC2'], prop=fontP

plt.ylabel('Depth (km)')
plt.gca().invert_yaxis()
plt.xscale('log')
plt.yticks(np.arange(0,13)/2.)

plt.savefig( 'boxplot_TMsinksensitivity_withavgdist.pdf', bbox_inches="tight")
plt.show()
