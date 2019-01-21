# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 15:31:22 2017

@author: nooteboom
"""

from parcels import (FieldSet, Field, ParticleSet, JITParticle,
                     ErrorCode, ParticleFile, Variable)
from datetime import timedelta as delta
from datetime import date
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
    """
    curr = start
    result = []
    result.append('{:%Y%m%d}'.format(curr))
    while curr <= end:
     #   yield curr
        curr += delta
        result.append('{:%Y%m%d}'.format(curr))
    del result[-1]
    result = np.array(result)
    return result

dd = 100. #The dwelling depth

dirwrite = '/OFES/OFESres/particledata/surface/'# Directory to write the data

posidx = int(sys.argv[1]) #ID of the file to define latitude and longitude ranges from 0 to 287

print 'id: ',posidx

#minpoints = np.mgrid[-180:179:30,-77:89:90]
#lonsmin = np.arange(0,359,20) - 180
#lonsmax = np.arange(20,361,20);#lonsmax[-1]=359;
#lonsmax = lonsmax - 180
#latsmin = np.array([-75,-50,-25,0, 25, 50])#np.array([-75,0])#np.array([-75,-25,0,25])
#latsmax = np.array([-50,-25,0,25, 50, 75])#np.array([-1,74])#np.array([-24,-1,24,74])
#Chaning the lon and lat you must also do this within the kernels
#minlon = lonsmin[posidx%18]
#maxlon = lonsmax[posidx%18]
#minlat = latsmin[posidx/18]
#maxlat = latsmax[posidx/18]


Bind = posidx

latsmin = np.arange(-75,75,3)
latsmax = np.arange(-72,76,3)

minlat = latsmin[posidx];maxlat = latsmax[posidx]
minlon=-180;maxlon=180;

latminind = (max(-75,minlat-1)+75)*10
latmaxind = (min(75,maxlat+1)+75)*10

res = 1 #resolution in degrees 

grid = np.mgrid[int(minlon):int(maxlon):res,int(minlat):int(maxlat):res]+0.5
n=grid[0].size
lons = np.reshape(grid[0],n)
lats = np.reshape(grid[1],n)


#Delete the particles on the land
landfilename = "grids/topography_OFES.nc"
Land = GetOFESLandArray(landfilename, 'HT')
[lonsz,latsz]=[np.array([lo for lo, la in zip(lons,lats) if Land[0,lo,la,0]>dd  ]),np.array([la for lo, la in zip(lons,lats) if Land[0,lo,la,0]>dd ])]

lonsz[lonsz<0] += 360

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

dep = dd * np.ones(latsz.shape) #2.5 m is minimum for OFES data

print 'lons: ',np.unique(lons)
print 'lats shape: ',lats.shape
print 'dwelling depth:', dd

sys.stdout.flush()
    
#%%
directory = '/projects/0/palaeo-parcels/OFESdata/'
def set_ofes_fieldset(snapshots):
    ufiles = [path.join(directory, 'y'+s[:4], 'u_vel', "nest_1_"+s+"000000u.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103
    vfiles = [path.join(directory, 'y'+s[:4], 'v_vel', "nest_1_"+s+"000000v.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103  
    tfiles = [path.join(directory, 'y'+s[:4], 'temp', "nest_1_"+s+"000000t.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103    

    sfiles = [path.join(directory, 'y'+s[:4], 'salinity', "nest_1_"+s+"000000s.nc".format(s)) for s in snapshots]
    

    filenames = { 'U': ufiles, 'V': vfiles,'temp': tfiles, 'salin': sfiles}
    variables = { 'U': 'zu', 'V': 'zv','temp': 'temperature', 'salin':'salinity'}
    
    dimensions = {      'U':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time'},
                        'V':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time'},
                        'temp':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time'},
                        'salin':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time'}}

    indices = {'lat': range(latminind, latmaxind)}

    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, indices=indices, allow_time_extrapolation=False)

    return fieldset

def SampleSurf(particle, fieldset, time, dt):
    particle.temp = fieldset.temp[time+dt, particle.lon, particle.lat,fieldset.surface]
    particle.salin = fieldset.salin[time+dt, particle.lon, particle.lat,fieldset.surface]             

def Age(particle, fieldset, time, dt):
    particle.age = particle.age + math.fabs(dt)  

def DeleteParticle(particle, fieldset, time, dt):
    particle.delete() 

def run_corefootprintparticles(dirwrite,outfile,lonss,lonsz,latss,latsz,dep):
    
    snapshots = snapshot_function(date(2000, 1, 3), date(2005, 12, 29), delta(days=3))
    fieldset = set_ofes_fieldset(snapshots)
    fieldset.add_periodic_halo(zonal=True)
    fieldset.add_constant('dwellingdepth', np.float(dd))
    fieldset.add_constant('maxage', 300000.*86400)
    fieldset.add_constant('surface', 2.6)

    class FixParticle(JITParticle):
        temp = Variable('temp', dtype=np.float32, initial=np.nan)
        age = Variable('age', dtype=np.float32, initial=0., to_write = False)
        salin = Variable('salin', dtype=np.float32, initial=np.nan)
  

    psets = ParticleSet.from_list(fieldset=fieldset, pclass=FixParticle, 
                                  lon=lonsz.tolist(), lat=latsz.tolist(),
                                     depth=dep, time = np.datetime64('2005-12-25'))#
    

    pfiles = ParticleFile(dirwrite + 'surface' + outfile, psets, outputdt=delta(days=3))

    kernelss =   psets.Kernel(SampleSurf) + psets.Kernel(Age) 

    psets.execute(kernelss, runtime=delta(days=2170), dt=delta(days=-1), output_file=pfiles, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})


outfile = "grid_dd" + str(int(dd)) + "_id"+str(posidx) +"_res"+str(res)
run_corefootprintparticles(dirwrite,outfile,lons,lonsz,lats,latsz,dep)

