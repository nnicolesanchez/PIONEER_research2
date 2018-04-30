# This script pulls out the following properties for the 
# ChaNGa galaxies: Patient 0, GM1, GM4, GM5, GM6
#        - Stellar mass of galaxy
#        - Gas mass of galaxy
#        * Split into CGM vs Disk?
#        - Merger? yes no
#            - Mass ratio if yes
#            - Redshift if yes

# N. Nicole Sanchez -- June 29 2017
# Univ. of Wash.    -- Nbody Shop
import matplotlib.pyplot as plt
import numpy as np
import pynbody
import sys

sims = ['/nobackupp8/fgoverna/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.004096','/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.004096','/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.004096','/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.004096','/nobackupp8/fgoverna/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.004096','/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.004096']
labels = ['P0','GM1','GM4','GM5','GM6','GM7']

import csv
with open('gxy_properties.csv', 'w') as csvfile:
    fieldnames = ['labels', 'total gas mass', 'total stellar mass', 'disk gas mass', 'disk stellar mass', 'CGM gas mass', 'CGM stellar mass', 'merger', 'merger mass ratio', 'merger redshift']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=',')

    writer.writeheader()

    for i in range(len(sims)):
        print(labels[i])
        sim = pynbody.load(sims[i],verbose=True)
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
        print('Total mass in gas in H1:',np.sum(h1.g['mass']))
        print('Total mass in stars in H1:',np.sum(h1.s['mass']))
    
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
    

        writer.writerow({'labels': str(labels[i]), 'total gas mass': str(np.sum(h1.g['mass'])), 'total stellar mass': str(np.sum(h1.s['mass'])), 'disk gas mass': str(np.sum(disk_gas['mass'])), 'disk stellar mass': str(np.sum(disk_stars['mass'])), 'CGM gas mass': str(np.sum(CGM_gas['mass'])), 'CGM stellar mass': str(np.sum(CGM_stars['mass'])), 'merger': 'Not found yet', 'merger mass ratio': 'Not found yet','merger redshift': 'Not found yet'})

        

