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

    return eta,etaz,vvir


prefix = '/home/christenc/Data/Sims/'
step = '00512'

dir799 = prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
file799 = 'h799.cosmo25cmb.3072g14HBWK'
key799 = 'h799'
dir516 = prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
file516 = 'h516.cosmo25cmb.3072g14HBWK'
key516 = 'h516'
dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
file986 = 'h986.cosmo50cmb.3072g14HBWK'
key986 = 'h986'
dir603 = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
file603 = 'h603.cosmo50cmb.3072g14HBWK'
key603 = 'h603'
dir258 = prefix + 'h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/'
file258 = 'h258.cosmo50cmb.3072g14HMbwK'
key258 = 'h258'
dir285 = prefix + 'h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/'
file285 = 'h285.cosmo50cmb.3072g14HMbwK'
key285 = 'h285'
dir239 = prefix + 'h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
file239 = 'h239.cosmo50cmb.3072g14HMbwK'
key239 = 'h239'
dirs    = np.array([ dir799, dir799, dir799, dir516, dir516, dir986, dir986, dir986, dir986, dir986, dir986, dir603, dir603, dir603, dir258, dir258, dir285, dir285, dir285, dir239])
files   = np.array([file799,file799,file799,file516,file516,file986,file986,file986,file986,file986,file986,file603,file603,file603,file258,file258,file285,file285,file285,file239])
haloid =  np.array([ 1   , 4   , 6   , 1   , 2   , 1   , 2   , 3   , 8   , 15  , 16  , 1   , 2    , 3  , 1   , 4   , 1   , 4   , 9   , 1   ])
masssort =np.array([10,      2,      9,      1,      8,      15,     18,     4,      0,      13,     17,     3,      7,      6,      12,     5,      11,     14,     16,     19])
dirs = dirs[masssort]
files = files[masssort]
haloid = haloid[masssort]

vvir = np.array(len(dirs))
eta_inner = np.array(len(dirs))
etaz_inner = np.array(len(dirs))
eta_outer = np.array(len(dirs))
etaz_outer = np.array(len(dirs))

for i in range(0,len(dirs)):
    filename = dirs[i] + files[i] + '.' + step + '/' + files[i] + '.' + step
    print(filename,haloid[i])
    radius = 0.25
    eta_inner,etaz_inner,vvir = calc_eta(filename,haloid[i],radius)
    print(eta_inner,etaz_inner)
    radius = 1
    eta_outer,etaz_outer,vvir = calc_eta(filename,haloid[i],radius)
    print(eta_outer,etaz_outer)    
