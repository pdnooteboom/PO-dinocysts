# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 15:31:22 2017

@author: nooteboom
"""

from parcels import (FieldSet, ParticleSet, JITParticle,
                     ErrorCode, ParticleFile, Variable)
from datetime import timedelta as delta
from datetime import date, datetime
import numpy as np
from os import path
import math
import sys
    
def snapshot_function(start, end, delta):
    """
    
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

directory = '/DATA_3D/snap_3day/'#folder with OFES flow field data
dirwrite = '/OFES/OFESres/wholetraj/'#Write directory for the output

#print maxlat, minlat, maxlon, minlon

lons = [154.91]#[120.5]#[147.79]#[154.91]#[148.8]
lats = [-49.71]#[23.5]#[-45.22]#[-49.26]

time = np.array([datetime(2005, 12, 25)]) 
lons = np.array(lons*time.shape[0])
lats = np.array(lats*time.shape[0])

print 'lons: ',np.unique(lons)
print 'lats shape: ',lats.shape

sys.stdout.flush()

#%%
def set_ofes_fieldset(snapshots):
    ufiles = [path.join(directory, 'y'+s[:4], 'u_vel', "nest_1_"+s+"000000u.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103
    vfiles = [path.join(directory, 'y'+s[:4], 'v_vel', "nest_1_"+s+"000000v.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103
    wfiles = [path.join(directory, 'y'+s[:4], 'w_vel', "nest_1_"+s+"000000w.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103   
    tfiles = [path.join(directory, 'y'+s[:4], 'temp', "nest_1_"+s+"000000t.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103    

    sfiles = [path.join(directory, 'y'+s[:4], 'salinity', "nest_1_"+s+"000000s.nc".format(s)) for s in snapshots]
    
    bfile = '/home/students/3830241/PhD/OFES_parcels/' + 'topography_OFES.nc'

    filenames = {'U': ufiles, 'V': vfiles, 'W': wfiles, 'temp': tfiles, 'salin': sfiles, 'B':bfile}
    variables = {'U': 'zu', 'V': 'zv', 'W': 'zw', 'temp': 'temperature', 'salin':'salinity', 'B': 'HT'}
    
    dimensions = {      'U':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'V':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'W':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'temp':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'salin':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'B':{'lat': 'LAT', 'lon': 'LONN1799_1800', 'time': 'TIME', 'depth': 'LEV'}} 
        
    indices = {'lat': range(0,700)}                 

    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, indices = indices, allow_time_extrapolation=False)

    return fieldset
        

def periodicBC(particle, fieldSet, time, dt):
    if particle.lon > 360:
        particle.lon -= 360        
    if particle.lon < 0:
        particle.lon += 360   
        
#Sink Kernel if only atsf is saved:
def Sink(particle, fieldset, time, dt):
    particle.temp = fieldset.temp[time+dt, particle.lon, particle.lat, fieldset.surface]
    particle.salin = fieldset.salin[time+dt, particle.lon, particle.lat, fieldset.surface]                     

def Age(particle, fieldset, time, dt):
    particle.age = particle.age + math.fabs(dt)  

def DeleteParticle(particle, fieldset, time, dt):
    particle.delete()

def initials(particle, fieldset, time, dt):
    if particle.age==0.:
        particle.depth = fieldset.surface
        particle.lon0 = particle.lon
        particle.lat0 = particle.lat
        particle.depth0 = particle.depth

def run_corefootprintparticles(dirwrite,outfile,lonss,latss,time):
    
    snapshots = snapshot_function(date(1980, 1, 3), date(2005, 12, 29), delta(days=3))
    fieldset = set_ofes_fieldset(snapshots)
    fieldset.add_periodic_halo(zonal=True)
    fieldset.B.allow_time_extrapolation = True
    fieldset.add_constant('maxage', 300000.*86400)
    fieldset.add_constant('surface', 2.6)

    class DinoOfesParticle(JITParticle):
        temp = Variable('temp', dtype=np.float32, initial=np.nan)
        age = Variable('age', dtype=np.float32, initial=0.)
        salin = Variable('salin', dtype=np.float32, initial=np.nan)
        lon0 = Variable('lon0', dtype=np.float32, initial=0., to_write=False)
        lat0 = Variable('lat0', dtype=np.float32, initial=0., to_write=False)
        depth0 = Variable('depth0',dtype=np.float32, initial=0., to_write=False)
        time0 = Variable('time0',dtype=np.float32, initial=0., to_write=False)
  

    pset = ParticleSet.from_list(fieldset=fieldset, pclass=DinoOfesParticle, lon=lonss.tolist(), lat=latss.tolist(), time = time) #depth=depths, 

    pfile = ParticleFile(dirwrite + outfile, pset, outputdt=delta(days=6))

    kernels = pset.Kernel(initials) + Sink + Age + periodicBC 

    pset.execute(kernels, runtime=delta(days=7300), dt=delta(days=-3), output_file=pfile, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
#days = 2555 for the 7 years run

outfile = "surfaceloc" + '_lon'+str(lons[0]) + '_lat' + str(lats[0])
run_corefootprintparticles(dirwrite,outfile,lons,lats,time)

print 'Exection finished'
