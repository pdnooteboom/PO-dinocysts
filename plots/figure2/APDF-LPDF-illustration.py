# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 11:38:49 2019

@author: nooteboom
"""

import numpy as np
from matplotlib.pylab import plt
import seaborn as sns
import matplotlib
from matplotlib.font_manager import FontProperties

hfont = {'fontname':'Helvetica'}

sns.set(style='darkgrid', font='Helvetica', font_scale=1.2)

fontP = FontProperties()
fontP.set_size('large')

sigma = 0.01

LPDF = np.random.normal(10,sigma,3000)
APDF1 = np.random.normal(15,sigma,3000)
APDF2 = np.append(np.random.normal(10,sigma,2000), np.random.normal(15,sigma,1000))
APDF3 = np.append(np.random.normal(10,sigma,2000), np.random.normal(5,sigma,1000))

#%%
fig = plt.figure(figsize=(13,3))
ax = fig.add_subplot(122)
sns.kdeplot(LPDF, bw=1.5, linewidth = 5, label='LPDF', color='navy');
sns.kdeplot(APDF1, bw=1.5, linewidth = 3, label='APDF1', color='firebrick', linestyle='--');
sns.kdeplot(APDF2, bw=1.5, linewidth = 3, label='APDF2', color='mediumseagreen', linestyle='--');
sns.kdeplot(APDF3, bw=1.5, linewidth = 3, label='APDF3', color='goldenrod', linestyle='--');
plt.xlabel('Temperature ($^{\circ}$C)', fontsize=18, **hfont)
plt.xlim(1,20)
for tick in ax.xaxis.get_major_ticks(): tick.label.set_fontsize(16) ;
for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(16) ;
plt.title('(b)', fontsize=18, **hfont)
#%%
xL = 0; yL = 0;
sigma = 6
bias = 20

res = 3

x1 = np.random.normal(xL+bias,sigma,3000)[::res]
x2 = np.append(np.random.normal(xL,sigma,2000) , np.random.normal(xL+bias,sigma,1000))[::res]
x3 = np.append(np.random.normal(xL,sigma,2000) , np.random.normal(xL-bias,sigma,1000))[::res]
y1 = np.random.normal(yL+bias,sigma,3000)[::res]
y2 = np.append(np.random.normal(yL,sigma,2000) , np.random.normal(yL+bias,sigma,1000))[::res]
y3 = np.append(np.random.normal(yL,sigma,2000) , np.random.normal(yL-bias,sigma,1000))[::res]

def temp(x,y):
    res = np.zeros(x.shape)
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            if(np.sqrt(x[i,j]**2 + y[i,j]**2) != 0):
                res[i,j] = 10+0.5*bias/80.*(x[i,j] + y[i,j])
                #/np.float(np.sqrt(x[i,j]**2 + y[i,j]**2))
            else: 
                res[i,j] = 10
    return res
lon, lat = np.mgrid[-40:40,-40:40]
T = temp(lon, lat)

sns.set_style("dark")
ax = fig.add_subplot(121)
plt.scatter(x1,y1,color='firebrick',alpha=0.15)
plt.scatter(x2,y2,color='mediumseagreen',alpha=0.15)
plt.scatter(x3,y3,color='goldenrod',alpha=0.15)
plt.scatter(xL,yL, color='navy')
levels = [2.5,7.5,12.5, 17.5]
CS = ax.contour(lon, lat, T, levels = [2.5,7.5,12.5, 17.5], linewidths=1, colors='k')
fmt = '%r $^{\circ}$C'
ax.clabel(CS, inline=1, fontsize=10, fmt=fmt, **hfont)
plt.xlabel('$^{\circ}$E', fontsize=18, **hfont)
plt.ylabel('$^{\circ}$N', fontsize=18, **hfont)
ax.tick_params(labelbottom=False, labelleft=False) 
plt.title('(a)', fontsize=18, **hfont)
plt.xlim(-40,40)
plt.ylim(-40,40)

plt.savefig('/Users/nooteboom/Documents/PhD/firstpaper/articleplots/plots/' + 'illustration_PDFS.pdf', bbox_inches="tight")
plt.show()