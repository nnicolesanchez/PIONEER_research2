# This script is a modified version of:

#Charlotte Christensen
#5/10/18
#This program will generate the mass loading and metal mass loading factors for a set of halos using the radial flux definition from Muratov 2015 and 2017

#Run with
#%run /home/christensen/Code/python/python_analysis/eta_metal_flux.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/eta_metal_flux.py
#ipython --pylab


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pynbody, sys
from pynbody.analysis import profile, angmom, halo
from pynbody import units, config
import pynbody.filt as f
import math
import os

def calc_eta(tfile,halo_num,radius):
    drad = 0.05
    zsolar = 0.0130215
    s = pynbody.load(filename)
    h = s.halos()
    if radius < (1 - drad):
        halo = h.load_copy(halo_num)
    else: 
        halo = h[halo_num]    
    halo.physical_units()

    pynbody.analysis.halo.center(halo)
    rvir = pynbody.array.SimArray(np.sqrt(np.max(halo['x'].in_units('kpc')**2 + halo['y'].in_units('kpc')**2 + halo['z'].in_units('kpc')**2)),'kpc')
    mvir = sum(halo['mass'].in_units('kg'))
    #vvir = ((units.G*mvir/rvir)**(1,2)).in_units('km s**-1')
    vvir = mvir #.in_units('kg')
    vvir = vvir/rvir.in_units('m')
    vvir = np.sqrt(6.67408e-11)*vvir**(0.5)/1000
    print(vvir)#in km/s
    radius1 = (radius  - drad)*rvir
    radius2 = (radius + drad)*rvir
        
    if radius < (1 - drad):
        annuli = halo[f.Annulus(radius1, radius2, cen=(0, 0, 0))].gas
    else:
        annuli = s[f.Annulus(radius1, radius2, cen=(0, 0, 0))].gas
    annuli['vrad'] = (annuli['x']*annuli['vx'] + annuli['y']*annuli['vy'] + annuli['z']*annuli['vz'])/np.sqrt(annuli['x']*annuli['x'] + annuli['y']*annuli['y'] + annuli['z']*annuli['z'])
    inward = pynbody.filt.LowPass('vrad','0 km s**-1')
    outward = pynbody.filt.HighPass('vrad','0 km s**-1')
    #pynbody.plot.sph.velocity_image(annuli[inward],vector_color="cyan",width = rvir*2,vector_resolution=100)
    #pynbody.plot.sph.velocity_image(annuli[outward],vector_color="cyan",width = rvir*2,vector_resolution=100)

    annuli['net_flux_mass_part'] =  (annuli['x']*annuli['vx'].in_units('kpc yr**-1') + annuli['y']*annuli['vy'].in_units('kpc yr**-1') + annuli['z']*annuli['vz'].in_units('kpc yr**-1'))/np.sqrt(annuli['x']*annuli['x'] + annuli['y']*annuli['y'] + annuli['z']*annuli['z'])*annuli['mass']/(radius2 - radius1)
    annuli['net_flux_z_part'] =  (annuli['x']*annuli['vx'].in_units('kpc yr**-1') + annuli['y']*annuli['vy'].in_units('kpc yr**-1') + annuli['z']*annuli['vz'].in_units('kpc yr**-1'))/np.sqrt(annuli['x']*annuli['x'] + annuli['y']*annuli['y'] + annuli['z']*annuli['z'])*annuli['mass']*(annuli['OxMassFrac']*2.09 + annuli['FeMassFrac']*1.06)/(radius2 - radius1)
    net_flux_mass = sum(annuli['net_flux_mass_part'])
    net_flux_z = sum(annuli['net_flux_z_part'])
    
    annuli[inward]['flux_mass_part'] =  (annuli[inward]['x']*annuli[inward]['vx'].in_units('kpc yr**-1') + annuli[inward]['y']*annuli[inward]['vy'].in_units('kpc yr**-1') + annuli[inward]['z']*annuli[inward]['vz'].in_units('kpc yr**-1'))/np.sqrt(annuli[inward]['x']*annuli[inward]['x'] + annuli[inward]['y']*annuli[inward]['y'] + annuli[inward]['z']*annuli[inward]['z'])*annuli[inward]['mass']/(radius2 - radius1)
    annuli[inward]['flux_z_part'] =  annuli[inward]['flux_mass_part']*(annuli[inward]['OxMassFrac']*2.09 + annuli[inward]['FeMassFrac']*1.06)
    annuli[outward]['flux_mass_part'] =  (annuli[outward]['x']*annuli[outward]['vx'].in_units('kpc yr**-1') + annuli[outward]['y']*annuli[outward]['vy'].in_units('kpc yr**-1') + annuli[outward]['z']*annuli[outward]['vz'].in_units('kpc yr**-1'))/np.sqrt(annuli[outward]['x']*annuli[outward]['x'] + annuli[outward]['y']*annuli[outward]['y'] + annuli[outward]['z']*annuli[outward]['z'])*annuli[outward]['mass']/(radius2 - radius1)
    annuli[outward]['flux_z_part'] =  annuli[outward]['flux_mass_part']*(annuli[outward]['OxMassFrac']*2.09 + annuli[outward]['FeMassFrac']*1.06)
    in_flux_mass = sum(annuli[inward]['flux_mass_part'])
    in_flux_z = sum(annuli[inward]['flux_z_part'])
    out_flux_mass = sum(annuli[outward]['flux_mass_part'])
    out_flux_z = sum(annuli[outward]['flux_z_part'])
    
    dtime = 50 #Myr 
    young = pynbody.filt.LowPass('age', str(dtime) + ' Myr')
    sfr = sum(halo.star[young]['mass'])/dtime/1e6
    
    eta = out_flux_mass/sfr
    etaz = out_flux_z/sfr

    return eta,etaz,vvir,in_flux_z,out_flux_z,halo.properties['time'].in_units('Gyr'),halo.properties['z']


#prefix = '/home/christenc/Data/Sims/'
#step = '00512'
prefix = '/nobackupp2/nnsanche/'
prefixnoBH = '/nobackupp2/nnsanche/NO_BHs/'

dirP0   = prefix + 'pioneer50h243.1536g1bwK1BH'
fileP0  = 'pioneer50h243.1536gst1bwK1BH'
keyP0   = 'P0'
dirGM1  = prefix + 'pioneer50h243GM1.1536gst1bwK1BH'
fileGM1 = 'pioneer50h243GM1.1536gst1bwK1BH'
keyGM1  = 'GM1'
dirGM7  = prefix + 'pioneer50h243GM7.1536gst1bwK1BH'
fileGM7 = 'pioneer50h243GM7.1536gst1bwK1BH'
keyGM7  = 'GM7'
dirGM4  = prefix + 'pioneer50h243GM4.1536gst1bwK1BH'
fileGM4 ='pioneer50h243GM4.1536gst1bwK1BH'
keyGM4  = 'GM4'
dirP0noBH   = prefixnoBH + 'pioneer50h243.1536gst1bwK1'
fileP0noBH  = 'pioneer50h243.1536gst1bwK1'
keyP0noBH   = 'P0noBH'
dirGM1noBH  = prefixnoBH + 'pioneer50h243GM1.1536gst1bwK1'
fileGM1noBH = 'pioneer50h243GM1.1536gst1bwK1'
keyGM1noBH  = 'GM1noBH'
dirGM7noBH  = prefixnoBH + 'pioneer50h243GM7.1536gst1bwK1'
fileGM7noBH = 'pioneer50h243GM7.1536gst1bwK1'
keyGM7noBH  = 'GM7noBH'
dirGM4noBH  = prefixnoBH + 'pioneer50h243GM4.1536gst1bwK1'
fileGM4noBH ='pioneer50h243GM4.1536gst1bwK1'
keyGM4noBH  = 'GM4noBH'

dirs  = np.array([dirP0, dirGM1, dirGM7, dirGM4, dirP0noBH, dirGM1noBH, dirGM7noBH, dirGM4noBH])
files = np.array([fileP0, fileGM1, fileGM7, fileGM4, fileP0noBH, fileGM1noBH, fileGM7noBH, fileGM4noBH])
haloid   = np.array([1,1,1,1,1,1,1,1])
labels = ['P0','GM1','GM7','GM4','P0noBH','GM1noBH','GM7noBH','GM4noBH']
#masssort = np.array([1,2,3,4,5,6,7,8])
#dirs = dirs[masssort]
#files = files[masssort]
#haloid = haloid[masssort]

vvir = np.array(len(dirs))
eta_inner = np.array(len(dirs))
etaz_inner = np.array(len(dirs))
eta_outer = np.array(len(dirs))
etaz_outer = np.array(len(dirs))

for i in range(3,len(dirs)):
    inner_influx_metal  = []
    inner_outflux_metal = []
    outer_influx_metal  = []
    outer_outflux_metal = []
    time     = []
    redshift = []
    if (os.path.exists(labels[i]+'_times.txt') == False):
        steps = np.loadtxt('../'+labels[i]+'/timesteps.txt',dtype=str)
        for ts in range(1,len(steps)-2):
            filename = dirs[i] + '/' + files[i] + '.00' + steps[ts]
            print(files[i],haloid[i],steps[ts])
            radius = 0.1
            eta_inner,etaz_inner,vvir,inner_influx_z,inner_outflux_z,t,red = calc_eta(filename,haloid[i],radius)
            radius = 1
            eta_outer,etaz_outer,vvir,outer_influx_z,outer_outflux_z,t,red = calc_eta(filename,haloid[i],radius)
            inner_influx_metal.append(inner_influx_z)
            inner_outflux_metal.append(inner_outflux_z)
            outer_influx_metal.append(outer_influx_z)
            outer_outflux_metal.append(outer_outflux_z)
            time.append(t)
            redshift.append(red)
    
            np.savetxt(labels[i]+'_inner_influx_metal.txt',inner_influx_metal)
            np.savetxt(labels[i]+'_inner_outflux_metal.txt',inner_outflux_metal)
            np.savetxt(labels[i]+'_outer_influx_metal.txt',outer_influx_metal)
            np.savetxt(labels[i]+'_outer_outflux_metal.txt',outer_outflux_metal)
            np.savetxt(labels[i]+'_times.txt',time)
            np.savetxt(labels[i]+'_redshifts.txt',redshift)

    else:
        inner_influx_metal = np.loadtxt(labels[i]+'_inner_influx_metal.txt')    
        inner_outflux_metal = np.loadtxt(labels[i]+'_inner_outflux_metal.txt')
        outer_influx_metal = np.loadtxt(labels[i]+'_outer_influx_metal.txt')    
        outer_outflux_metal = np.loadtxt(labels[i]+'_outer_outflux_metal.txt')
        time = np.loadtxt(labels[i]+'_times.txt')
        redshift = np.loadtxt(labels[i]+'_redshifts.txt')

    plt.plot(time,outer_influx_metal,label='Inflow at Rvir',linestyle='-',linewidth=2,color='SteelBlue')
    plt.plot(time,outer_outflux_metal,label='Outflow at Rvir',linestyle=':',linewidth=2,color='SteelBlue')    
    plt.plot(time,inner_influx_metal,label='Inflow at 0.1*Rvir',linestyle='-',linewidth=2,color='SkyBlue')
    plt.plot(time,inner_outflux_metal,label='Outflow at 0.1*Rvir',linestyle=':',linewidth=2,color='SkyBlue')
    plt.ylim(-1,1)
    plt.xlim(2,14)
    plt.title(labels[i])
    plt.legend()
    plt.savefig(labels[i]+'_in_outflow_metalmass_time.pdf')
#    plt.show()
    plt.clf()
