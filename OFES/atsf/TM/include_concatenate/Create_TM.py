#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Script to compute Transition matrices from given simulated trajectories

Note:
    - full indices' refer to the whole area including land, 'reduced indnices' without land

@author: wichmann
"""

from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import os


#pltse = 'annual' 

#dt = -1 #Time for loading particles (TM timestep)
#ddeg=3. #Spacing of TM

#Lons = np.arange(280,325.,ddeg)
#Lats = np.arange(25.,47.,ddeg)


class Particle(object):
    def __init__(self,lon,lat, Lons, Lats):
        self.lon=lon
        self.lat=lat
        self.surv= [False if (np.isnan(lon[j]) or lon[j]>max(Lons) or lon[j]<min(Lons) or lat[j]>max(Lats) or lat[j]<min(Lats)) else True for j in range(len(lat))] #Label for non-deleted particles#or lon[j]>max(Lons) or lon[j]<min(Lons) or


class TM_ParticleSet(object):
    
    def __init__(self,lon,lat, Lons, Lats):
        print 'Total number of particles: ', len(lon)
        self.all_particles = np.empty(len(lon),dtype = Particle) #Container for all particiles (including those that aare killed)
        self.survived =[]
        for i in range(len(self.all_particles)):
            if (i%20000==0):
                print 'prop. set up particles: ', i/np.float(len(self.all_particles))
            self.all_particles[i]=Particle(lon[i],lat[i], Lons, Lats)
            self.survived.append(self.all_particles[i].surv)
        
        self.lons_I=[lon[i][0] for i in range(len(lon))] #For computing initial particle number per cell later (we do not only consider those that survive the entire period)
        self.lats_I=[lat[i][0] for i in range(len(lon))]
        self.Lons = Lons
        self.Lats = Lats
    
    def setup_TM(self,ddeg,Lons,Lats):
        surv = [self.survived[i][1] for i in range(len(self.all_particles))]
        self.particles = self.all_particles[surv] #Take only particles that are not deleted

        self.Lons = Lons#np.arange(280,360.,ddeg)
        self.Lats = Lats #np.arange(0.,60.,ddeg)
        self.N_initial_I=np.zeros(self.Lons.size*self.Lats.size)         
        
        #Initial and final points for non-deleted particles
        lons_i = np.array([p.lon[0] for p in self.particles])
        lats_i = np.array([p.lat[0] for p in self.particles])
        lons_f = np.array([p.lon[1] for p in self.particles])
        lats_f = np.array([p.lat[1] for p in self.particles])
        print 'Number of non-deleted particles: ', len(lons_f)

        print 'Computing TM indices'
        self.TMindex_i=np.array([int(((la-np.min(self.Lats))//ddeg)*len(self.Lons)+(lo-np.min(self.Lons))//ddeg) for la,lo in zip(lats_i,lons_i)])
        self.TMindex_f=np.array([int(((la-np.min(self.Lats))//ddeg)*len(self.Lons)+(lo-np.min(self.Lons))//ddeg) for la,lo in zip(lats_f,lons_f)])
        self.TMindex_I=np.array([int(((la-np.min(self.Lats))//ddeg)*len(self.Lons)+(lo-np.min(self.Lons))//ddeg) for la,lo in zip(self.lats_I,self.lons_I)])
        
        print 'max i: ', self.TMindex_i[self.TMindex_i>4800]
        print 'max f: ', self.TMindex_f[self.TMindex_f>4800]
        
        
        print 'lon: ', lons_i[lons_i>325.]
        print 'lon: ', lons_i[lons_i<280.]

        print 'lat: ', lats_i[lats_i>47.]
        print 'lat: ', lats_i[lats_i<25.]

        
        print 'Nshape: ', self.N_initial_I.shape
        print 'Computing initial position'
        for i in range(len(self.lons_I)):
            if not np.isnan(self.lons_I[i]) or ma.is_masked(self.lons_I[i]):            
                self.N_initial_I[self.TMindex_I[i]]+=1.       
   
        #Identify the cells where not the max number of particles was put initially
        n=np.max(self.N_initial_I)
        self.coast = [1 if (self.N_initial_I[i]<n and self.N_initial_I[i]>0.) else 0 for i in range(len(self.N_initial_I))]
    
        print 'Compute TM'
        self.TM = np.zeros((len(self.N_initial_I),len(self.N_initial_I)))
        for i in range(len(self.particles)):
            i_start = self.TMindex_i[i]
            i_finish = self.TMindex_f[i]
            self.TM[i_finish,i_start]+=1./self.N_initial_I[self.TMindex_i[i]]
    
    def save_TM(self,name):
        np.savez(name, TM=self.TM, Lons=self.Lons, Lats=self.Lats, Ninit = self.N_initial_I, coast=self.coast)


def TM(datadir,filenames,outputdir):
    #Create TM for a specific time
    
    i = 0
    pfile = datadir +filenames + str(i) + '.nc'
    data = Dataset(pfile, 'r')
    
    lon=np.concatenate((data.variables['lons0'][:],data.variables['lon'][:]),axis=1)
    lat=np.concatenate((data.variables['lats0'][:],data.variables['lat'][:]),axis=1)#data.variables['lat'][:,[0,dt]]

    
    P=TM_ParticleSet(lon,lat)
    P.setup_TM(ddeg)
    P.save_TM(outputdir + 'TM_bin'+str(ddeg))


#def get_annual_TM(TMtotalName,TMfiles):
#    #Create annual TM from single time TMs
#    for t in range(8):
#        print 'Time: ', t
#        data=np.load(TMfiles + '_time'+str(t)+'.npz')
#        if t==0:
#            TM = data['TM']
#            Lons = data['Lons']    
#            Lats = data['Lats']
#            Coasts = data['coast']
#        else:
#            tm = data['TM']
#            TM=np.dot(tm,TM)
#    
#    np.savez(TMtotalName, TM=TM, Lons=Lons, Lats=Lats, coast=Coasts)



