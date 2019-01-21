

"""Calculate for every gridbox the mean (temperature) of the advected 
particles of every grid box in TM and the mean (temperature) at the fixed surface
location of the grid box """
from __future__ import division
import numpy as np
from netCDF4 import Dataset
import time
from numba import jit

#%% Some functions
from datetime import date, datetime, timedelta

Y = 2000 # dummy leap year to allow input X-02-29 (leap day)
seasons = [('winter', (date(Y,  1,  1),  date(Y,  3, 20))),
           ('spring', (date(Y,  3, 21),  date(Y,  6, 20))),
           ('summer', (date(Y,  6, 21),  date(Y,  9, 22))),
           ('autumn', (date(Y,  9, 23),  date(Y, 12, 20))),
           ('winter', (date(Y, 12, 21),  date(Y, 12, 31)))]

def get_season(now):
    if isinstance(now, datetime):
        now = now.date()
    now = now.replace(year=Y)
    return next(season for season, (start, end) in seasons
                if start <= now <= end)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
    
def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx    
    
def find_down(array,value):
    if(value>=array[0]):
        array = [n-value for n in array]
        idx = np.array([n for n in array if n<0]).argmax()
    else:
        idx = np.nan
    return idx     

#%% To open a file
res = 1
sp = 3.
print 'sp: ',sp
dd = 10. 

var = 'salin'
print 'var: ',var
ddeg = 1

dirRead = '/OFES/OFESres/particledata/'
dirWrite = '/OFES/OFESres/TM/box-meandiff/'

#%% Getting the average and variance differences 
 
pltse = 'annual' 
if(pltse=='winter'):
    se = 0
elif(pltse=='spring'):
    se = 1
elif(pltse=='summer'):
    se = 2
elif(pltse=='autumn'):
    se = 3
else:
    se = 4

lons0 = np.array([])
lats0 = np.array([])
lons = np.array([])
lats = np.array([])
season = np.array([])
temp = np.array([])

#Arrays for fixed surface locations
fixlon = np.array([])
fixlat = np.array([])
fixtemp = np.array([])

start = time.time()
for posidx in range(50):

    if(posidx==84 or posidx==85):
        print'no posidx file (because only land in this area):'
        print posidx
    else:
        if(posidx%10==0):
            print posidx
        nc = Dataset(dirRead + 'sp%d_dd%d/grid_id'%(int(sp),int(dd)) +str(posidx)+'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res) + ".nc") 
        ncn = Dataset(dirRead + 'surface/surfacegrid' +'_dd'+str(int(dd))+'_id'+str(posidx)+"_res"+str(res) + '.nc')
        
        if(nc['lon0'][0,0]<0):
            lons0 = np.append(lons0,nc['lon0'][:,0]+360)   
        else:
            lons0 = np.append(lons0,nc['lon0'][:,0])         
#        lons0 = np.concatenate((lons0,nc['lon0'][:,0]),axis=0)
        lats0 = np.concatenate((lats0,nc['lat0'][:,0]), axis=0)
        lons = np.concatenate((lons,nc['lon'][:,0]),axis=0)
        lats = np.concatenate((lats,nc['lat'][:,0]), axis=0) 
        temp = np.concatenate((temp,nc[var][:,0]), axis=0) 
        
        fixlon = np.concatenate((fixlon,ncn['lon'][:,-1]+360),axis=0)
        fixlat = np.concatenate((fixlat,ncn['lat'][:,-1]), axis=0) 
        fixtemp = np.concatenate((fixtemp,np.nanmean(ncn[var][:], axis=1)), axis=0)         

print 'time to load the data (minutes): ', ((time.time()-start)/60.)    

lons = lons%360
lons0 = lons0%360
fixlon = fixlon%360

print 'Are the longitudes in similar: '
print fixlon[0]
print lons[0]
print lons0[0]  

minlon = min(lons0)
maxlon = max(lons0)+1
minlat = min(lats0)
maxlat = max(lats0)+1
#%% Now loop over all regions to define the T matrix
    
#Define the grid at the bottom:
Lons = np.arange(minlon,maxlon+ddeg,ddeg)-0.5
Lats = np.arange(minlat,maxlat+ddeg,ddeg)-0.5
Lonss, Latss = np.meshgrid(Lons, Lats)

vLons = Lonss.reshape(Lonss.size) 
vLats = Latss.reshape(Latss.size)   

#print fixlon.shape
#print fixlat.shape
#print fixtemp.shape
#print fixlon
#print Lons
def region_index(lats0, lons0, fixlon, fixlat,  se):
    """Return the indices of the grid box at the bottom (Bi)""" 
    Bi = np.zeros(lons0.shape[0]) #The grid cell index where a particle is released at the bottom
    Si = np.zeros(fixlon.shape[0])
    
    #Look in which temperature region the particle surfaces and grid box region
    # released
    for i in range(lons0.shape[0]):

        if(pltse == 'annual'):#season[i]==se or se==4):
            lo = find_down(Lons,lons0[i]); lon = Lons[lo];
            la = find_down(Lats,lats0[i]); lat = Lats[la];
            
            Bi[i] = np.where(np.logical_and(vLons==lon,vLats==lat))[0]      

    for i in range(fixlon.shape[0]):

        if(pltse == 'annual'):#season[i]==se or se==4):
            lo = find_down(Lons,fixlon[i]); lon = Lons[lo];
            la = find_down(Lats,fixlat[i]); lat = Lats[la];
            
            Si[i] = np.where(np.logical_and(vLons==lon,vLats==lat))[0]  
            
    return Bi, Si
    
@jit(nopython=True)
def region_index_jit(lats0, lons0, fixlon, fixlat,  vLons, vLats, Lons, Lats):
    """Return the indices of the grid box at the bottom (Bi)""" 
    le = len(lons0)
    sle = len(fixlon)
    Bi = [0]*le #The grid cell index where a particle is released at the bottom
    Si = [0]*sle
    
    #Look in which temperature region the particle surfaces and grid box region
    # released
    for i in range(le):
        Lonsa = [n-lons0[i] for n in Lons]
        Lonsa = [n for n in Lonsa if n<0]
        ma = -10000.
        for j in range(len(Lonsa)):
            if(Lonsa[j]>ma):
                ma = Lonsa[j]
                lo = j
        lon = Lons[lo];
        
        Latsa = [n-lats0[i] for n in Lats]
        Latsa = [n for n in Latsa if n<0]
        ma = -10000.
        for j in range(len(Latsa)):
            if(Latsa[j]>ma):
                ma = Latsa[j]
                la = j
        lat = Lats[la];
           
        for j in range(len(vLons)):
            if(vLons[j]==lon):
                if(vLats[j]==lat):
                    Bi[i] = j   

    for i in range(sle):
        Lonsa = [n-fixlon[i] for n in Lons]
        Lonsa = [n for n in Lonsa if n<0]
        ma = -10000.
        for j in range(len(Lonsa)):
            if(Lonsa[j]>ma):
                ma = Lonsa[j]
                lo = j
        lon = Lons[lo];
        
        Latsa = [n-fixlat[i] for n in Lats]
        Latsa = [n for n in Latsa if n<0]
        ma = -10000.
        for j in range(len(Latsa)):
            if(Latsa[j]>ma):
                ma = Latsa[j]
                la = j
        lat = Lats[la];
           
        for j in range(len(vLons)):
            if(vLons[j]==lon):
                if(vLats[j]==lat):
                    Si[i] = j                      
                    
    return Bi, Si

#start = time.time()   
#Bi, Si =  region_index(lats0, lons0, fixlon, fixlat,  se)
#print 'time \'region_index\' (minutes): ', ((time.time()-start)/60.)

start = time.time()   
Bi, Si =  region_index_jit(lats0, lons0, fixlon, fixlat,  vLons, vLats, Lons, Lats)
print 'time \'region_index\' (minutes): ', ((time.time()-start)/60.)

#print 'Bi the same: ', (Bi==Bi2).all()
#print 'Si the same: ', (Si==Si2).all()
# %%Calculating the transition matrix: 
print 'Calculate the averages of every grid box'  
trm = np.full((len(vLons)),np.nan) 

#Function which returns both the mean of the fixed grid boxes and 
#the mean (temperature) of the particles which reach the surface after advection
@jit(nopython=True)
def makeTM(Bi, Si, lons0, temp, fixtemp):
    trm = [0.] * len(vLons)
    trmfix = [0.] * len(vLons)
    for r1 in range(len(vLons)):
        ts = 0.
        tsl = 0.
        for j in range(lons0.shape[0]):
            if(r1==Bi[j]):
                ts += temp[j] 
                tsl += 1.
        if(tsl>0):
            trm[r1] = ts / float(tsl)

        ts = 0.
        tsl = 0.
        for j in range(fixlon.shape[0]):
            if(r1==Si[j]):
                ts += fixtemp[j] 
                tsl += 1.
        if(tsl>0):
            trmfix[r1] = ts / float(tsl)        
 
    return trm, trmfix

start = time.time()

trm, trmfix = makeTM(Bi, Si, lons0, temp, fixtemp)
print 'time \'rTM calculation\' (minutes): ', ((time.time()-start)/60.)

#print trm
#print trmfix

if (pltse == 'annual'):
    np.savez(dirWrite + 'TM_box-meandiff'+var+'_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd)), TM=np.array(trm), TMfix = np.array(trmfix), Lons=Lons, Lats=Lats) 
else:
    np.savez(dirWrite +  'TM_box-meandiff'+var+'_ddeg%d_sp%d_dd%d_'%(ddeg,int(sp),int(dd)) + 'season' + pltse, TM=np.array(trm), TMfix = np.array(trmfix), Lons=Lons, Lats=Lats) 
          
        
