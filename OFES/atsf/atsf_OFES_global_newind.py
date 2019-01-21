# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 15:31:22 2017

@author: nooteboom
"""

from parcels import (FieldSet, Field, ParticleSet, JITParticle, AdvectionRK4_3D,
                     ErrorCode, ParticleFile, Variable)
from datetime import timedelta as delta
from datetime import date, datetime
import numpy as np
from os import path
import math
from netCDF4 import Dataset
import sys

def GetOFESLandArray(filename, fieldname):
    """
    Function to return a Field with 1's at land points and 0's at ocean points, based on python basemap
    :param f: a field of the .nc file, but not a parcels field object! For OFES, land points are masked. This is used here!
    """
    pfile = Dataset(filename, 'r')
    Lon = pfile.variables['LONN1799_1800'][:]
    Lat = pfile.variables['LAT'][:]
    f = pfile.variables[fieldname][:]
    f = f[0,0,:,:]    
    print f.shape, Lon.shape, Lat.shape#, L.shape   
    Land=Field('Land',f,transpose=False, lon=Lon,lat=Lat)
    return Land
    
def snapshot_function(start, end, delta):
    """
    This function returns a list of all fields which have to be loaded between 
    'start' end 'end' with 'delta' as timestep.
    ''
    """
    curr = start
    result = []
    result.append('{:%Y%m%d}'.format(curr))
    while curr <= end:
        curr += delta
        result.append('{:%Y%m%d}'.format(curr))
    del result[-1]
    result = np.array(result)
    return result

sp = 6. #The sinkspeed m/day
dd = 10. #The dwelling depth

directory = '' #Data of OFES flow/temperature/salinity field
dirwrite = '/OFES/OFESres/particledata/sp%d_dd%d/'%(int(sp), int(dd))# directory to write the data

posidx = int(sys.argv[1]) #ID of the file to define latitude and longitude ranges from 0 to 287

print 'id: ',posidx

#Chaning the lon and lat you must also do this within the kernels
Bind = posidx

latsmin = np.arange(-75,75,3)
latsmax = np.arange(-72,76,3)

minlat = latsmin[posidx];maxlat = latsmax[posidx]
minlon=-180;maxlon=180;

latminind = (max(-75,minlat-8)+75)*10
latmaxind = (min(75,maxlat+8)+75)*10

res = 1 #resolution of particle release in degrees 

grid = np.mgrid[int(minlon):int(maxlon):res,int(minlat):int(maxlat):res]+0.5
n=grid[0].size
lons = np.reshape(grid[0],n)
lats = np.reshape(grid[1],n)

#Delete the particles on the land
landfilename =  "grids/topography_OFES.nc"
Land = GetOFESLandArray(landfilename, 'HT')
[lonsz,latsz]=[np.array([lo for lo, la in zip(lons,lats) if Land[0,lo,la,0]>dd  ]),np.array([la for lo, la in zip(lons,lats) if Land[0,lo,la,0]>dd ])]

lonsz[lonsz<0] = lonsz[lonsz<0]+360

if(not lonsz.size):
    sys.exit("Only land in the run with this idx")
    
#Plot particles and density to check if it is correct
#from mpl_toolkits.basemap import Basemap
#import matplotlib.pyplot as plt
#fig = plt.figure(figsize=(22, 16))
#ax = fig.add_subplot(211)
#ax.set_title("Particles")
#m = Basemap(projection='merc', llcrnrlat=minlat,  urcrnrlat=maxlat, llcrnrlon=minlon, urcrnrlon=maxlon, resolution='l')
#m.drawparallels(np.array([20,40,60]), labels=[True, False, False, True])
#m.drawmeridians(np.array([280,300,320,340,360,380]), labels=[False, False, False, True])
#m.drawcoastlines()
#xs, ys = m(lons, lats)
#m.scatter(xs,ys)
#plt.show()



times = np.array([datetime(2005, 12, 25) - delta(days=x) for x in range(0,int(365*5+1),5)]) 
time = np.empty(shape=(0));lons = np.empty(shape=(0));lats = np.empty(shape=(0));
for i in range(len(times)):
    lons = np.append(lons,lonsz)
    lats = np.append(lats, latsz)
    time = np.append(time, np.full(len(lonsz),times[i])) 

dep = dd * np.ones(latsz.shape)

print 'lats: ',np.unique(lats)
print 'lats shape: ',lats.shape

sys.stdout.flush()

print 'sinking velocity: ', sp
print 'dwelling depth:', dd
print 'grid'

    
#%%
def set_ofes_fieldset(snapshots):
    ufiles = [path.join(directory, 'y'+s[:4], 'u_vel', "nest_1_"+s+"000000u.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103
    vfiles = [path.join(directory, 'y'+s[:4], 'v_vel', "nest_1_"+s+"000000v.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103
    wfiles = [path.join(directory, 'y'+s[:4], 'w_vel', "nest_1_"+s+"000000w.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103   
    tfiles = [path.join(directory, 'y'+s[:4], 'temp', "nest_1_"+s+"000000t.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103    

    sfiles = [path.join(directory, 'y'+s[:4], 'salinity', "nest_1_"+s+"000000s.nc".format(s)) for s in snapshots]
    
    bfile = 'grids/topography_OFES.nc'

    filenames = {'U': ufiles, 'V': vfiles, 'W': wfiles, 'temp': tfiles, 'salin': sfiles, 'B':bfile}
    variables = {'U': 'zu', 'V': 'zv', 'W': 'zw', 'temp': 'temperature', 'salin':'salinity', 'B': 'HT'}
    
    dimensions = {      'U':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'V':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'W':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'temp':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'salin':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'B':{'lat': 'LAT', 'lon': 'LONN1799_1800', 'time': 'TIME', 'depth': 'LEV'}}
        
    indices = {'lat': range(latminind, latmaxind)}                 

    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, indices = indices, allow_time_extrapolation=False)

    return fieldset
        

def periodicBC(particle, fieldSet, time, dt):
    if particle.lon > 360:
        particle.lon -= 360        
    if particle.lon < 0:
        particle.lon += 360   
        
#Sink Kernel if only atsf is saved:
def Sink(particle, fieldset, time, dt):
    if(particle.depth>fieldset.dwellingdepth):
        particle.depth = particle.depth + fieldset.sinkspeed * dt
    if(particle.depth<fieldset.dwellingdepth):
        particle.depth = fieldset.surface
        particle.temp = fieldset.temp[time+dt, particle.lon, particle.lat, particle.depth]
        particle.salin = fieldset.salin[time+dt, particle.lon, particle.lat, particle.depth]
        particle.delete()        

def SampleSurf(particle, fieldset, time, dt):
    particle.temp = fieldset.temp[time+dt, particle.lon, particle.lat, fieldset.surface]
    particle.salin = fieldset.salin[time+dt, particle.lon, particle.lat, fieldset.surface]              

def Age(particle, fieldset, time, dt):
    particle.age = particle.age + math.fabs(dt)  

def DeleteParticle(particle, fieldset, time, dt):
    particle.delete()

def initials(particle, fieldset, time, dt):
    if particle.age==0.:
        particle.depth = fieldset.B[time+dt, particle.lon, particle.lat, particle.depth]
        particle.lon0 = particle.lon
        particle.lat0 = particle.lat
        particle.depth0 = particle.depth
        
def FirstParticle(particle, fieldset, time, dt):
    particle.lon = particle.lon0
    particle.lat = particle.lat0
    particle.depth = fieldset.dwellingdepth   

def run_corefootprintparticles(dirwrite,outfile,lonss,lonsz,latss,latsz,dep,time): 
    snapshots = snapshot_function(date(2000, 1, 3), date(2005, 12, 29), delta(days=3))
    fieldset = set_ofes_fieldset(snapshots)
    fieldset.add_periodic_halo(zonal=True)
    fieldset.B.allow_time_extrapolation = True
    fieldset.add_constant('dwellingdepth', np.float(dd))
    fieldset.add_constant('sinkspeed', sp/86400.)
    fieldset.add_constant('maxage', 300000.*86400)
    fieldset.add_constant('surface', 2.6)

    class DinoOfesParticle(JITParticle):
        temp = Variable('temp', dtype=np.float32, initial=np.nan)
        age = Variable('age', dtype=np.float32, initial=0.)
        salin = Variable('salin', dtype=np.float32, initial=np.nan)
        lon0 = Variable('lon0', dtype=np.float32, initial=0.)
        lat0 = Variable('lat0', dtype=np.float32, initial=0.)
        depth0 = Variable('depth0',dtype=np.float32, initial=0., to_write=False)
        time0 = Variable('time0',dtype=np.float32, initial=0., to_write=False)

    pset = ParticleSet.from_list(fieldset=fieldset, pclass=DinoOfesParticle, lon=lonss.tolist(), lat=latss.tolist(), time = time) #depth=depths, 

    pfile = ParticleFile(dirwrite + outfile, pset, write_ondelete=True)

    kernels = pset.Kernel(initials) + Sink  + pset.Kernel(AdvectionRK4_3D) + Age + periodicBC 

    pset.execute(kernels, runtime=delta(days=2170), dt=delta(minutes=-60), output_file=pfile, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

outfile = "grid_id"+str(posidx) +'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res)
run_corefootprintparticles(dirwrite,outfile,lons,lonsz,lats,latsz,dep,time)

print 'Exection finished'
