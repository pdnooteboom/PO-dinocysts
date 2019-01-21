# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 20:56:46 2018

@author: nooteboom
"""

import numpy as np
import matplotlib.pylab as plt
import matplotlib
from matplotlib.font_manager import FontProperties
import seaborn as sns

sns.set(style="darkgrid")

font = {'family' : 'Helvetica',
#    'weight' : 'bold',
    'size'   : 22}

matplotlib.rc('font', **font)

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

#%%
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
plt.legend(bbox_to_anchor=(1.35, 1))#,['SC1', 'SC2'], prop=fontP

#plt.legend()
plt.xlabel('Sinking speed (m day$^{-1}$)')
plt.ylabel('Depth (km)')
plt.gca().invert_yaxis()
#plt.yscale('log')
plt.xscale('log')
#plt.xlim(0,70)
plt.yticks(np.arange(0,13)/2.)
#plt.xticks([3,6,11,25,50,100,200,500])
plt.savefig('/Users/nooteboom/Documents/PhD/firstpaper/articleplots/plots/' + 'illustration_sinkingspeeds.pdf', bbox_inches="tight")
plt.show()