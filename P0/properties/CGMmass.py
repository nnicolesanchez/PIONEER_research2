# This script find the mass of the CGM for a ChaNGa galaxy snapshot
#     - Outputs:
#          total stellar mass of disk across timesteps
#          total gas mass of disk across timesteps
#          total stellar mass of CGM across timesteps
#          total gas mass of CGM across timesteps
#      * Also prints out timesteps times in Gyr

# N. Nicole Sanchez -- June 29 2017
# Univ. of Wash.    -- Nbody Shop
import matplotlib.pyplot as plt
import numpy as np
import pynbody
import sys

#if len(sys.argv) == 1:
#    print('Syntax: python patient0_scalelength.py max_age min_age r_max z_max')
#    age_max = 6.
#else:
#    age_max = float(sys.argv[1])

timesteps = np.loadtxt('../timesteps.txt',dtype='str')

disk_star_mass = np.zeros(len(timesteps))
disk_gas_mass  = np.zeros(len(timesteps))
CGM_star_mass  = np.zeros(len(timesteps))
CGM_gas_mass   = np.zeros(len(timesteps))
redshifts      = np.zeros(len(timesteps))

#i = 21
for i in range(len(timesteps)):
    print('Starting timestep: ', timesteps[i])
    sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00'+timesteps[i],verbose=True)
    sim.properties
    sim.properties['time'].in_units('Gyr')
    sim.loadable_keys()
    sim.physical_units()
    h = sim.halos()
    h1 = h[1]
    pynbody.analysis.angmom.faceon(h1)

# Want to isolate CGM  
# Isolate and remove disk stars within radius 0-10 kpc & vertically 4 kpc 
    r_max = 10  # kpc
    z_max = 4   # kpc

# R = sqrt(x^2 + y^2)
    R_d = ((h1.s['x'].in_units('kpc'))**2. + (h1.s['y'].in_units('kpc'))**2.)**(0.5)
    disk_stars_xymax =  (R_d < r_max)
    disk_stars_zmax  = (h1.s['z'].in_units('kpc') < z_max) & (h1.s['z'].in_units('kpc') > -z_max)
    print('Removing H1 stars limited to within 10  kpc radially and 4 kpc vertically.')

    disk_stars_mask = disk_stars_xymax & disk_stars_zmax
    disk_stars = h1.s[disk_stars_xymax & disk_stars_zmax]
    CGM_stars  = h1.s[~disk_stars_mask]
    print('Total number of stars in H1: ', len(h1.s))
    print('Number of stars in disk: ', len(disk_stars))
    print('Mass of stars in disk: ', np.sum(disk_stars['mass']))
    print('Number of stars in CGM/Halo: ', len(CGM_stars))
    print('Mass of stars in CGM/Halo: ', np.sum(CGM_stars['mass']))


    Rg_d = ((h1.g['x'].in_units('kpc'))**2. + (h1.g['y'].in_units('kpc'))**2.)**(0.5)
    disk_gas_xymax =  (Rg_d < r_max)
    disk_gas_zmax  = (h1.g['z'].in_units('kpc') < z_max) & (h1.g['z'].in_units('kpc') > -z_max)

    disk_gas_mask = disk_gas_xymax & disk_gas_zmax
    disk_gas = h1.g[disk_gas_xymax & disk_gas_zmax]
    CGM_gas  = h1.g[~disk_gas_mask]
    print('Total number of gas particles in H1: ', len(h1.g))
    print('Number of gas particles in disk: ', len(disk_gas))
    print('Mass of gas in disk: ', np.sum(disk_gas['mass']))
    print('Number of gas particles in CGM/Halo: ', len(CGM_gas))
    print('Mass of gas in CGM/Halo: ', np.sum(CGM_gas['mass']))

    np.savetxt('mr_CGM_starmass_gasmass.txt',(np.sum(CGM_stars['mass']),np.sum(CGM_gas['mass'])))
    np.savetxt('mr_disk_starmass_gasmass.txt',(np.sum(disk_stars['mass']),np.sum(disk_gas['mass'])))
    
    redshifts[i]      = sim.properties['time'].in_units('Gyr')
    disk_star_mass[i] = np.sum(disk_stars['mass'])
    disk_gas_mass[i]  = np.sum(disk_gas['mass'])
    CGM_star_mass[i]  = np.sum(CGM_stars['mass'])
    CGM_gas_mass[i]   = np.sum(CGM_gas['mass'])

np.savetxt('disk_star_mass.txt',disk_star_mass)
np.savetxt('disk_gas_mass.txt',disk_gas_mass)
np.savetxt('CGM_star_mass.txt',CGM_star_mass)
np.savetxt('CGM_gas_mass.txt',CGM_gas_mass)
np.savetxt('redshifts.txt',redshifts)

pynbody.analysis.angmom.sideon(h1)
pynbody.plot.image(disk_stars,width='500 kpc', av_z=True, cmap='Blues')
plt.title('Stars in CGM of H1')
plt.savefig('h1_CGMstars_sideon.pdf')
plt.show()
plt.clf()

pynbody.analysis.angmom.faceon(h1)
pynbody.plot.image(CGM_stars,width='500 kpc', av_z=True, cmap='Blues')
plt.title('Stars in CGM of H1')
plt.savefig('h1_CGMstars_faceon.pdf')
plt.show()
plt.clf()

