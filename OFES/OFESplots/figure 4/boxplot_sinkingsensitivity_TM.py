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
plt.savefig( 'boxplot_TMsinksensitivity_withavgdist.pdf', bbox_inches="tight")
plt.show()
