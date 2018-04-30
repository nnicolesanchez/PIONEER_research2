import matplotlib.pyplot as plt
#import pandas as pd
import numpy as np
import pynbody

N = 10

GM1h1z0igasorder = np.loadtxt('../compare_GM1_GM4/GM1_h1_tracedigasorder.txt')
print(len(GM1h1z0igasorder))

timesteps = np.loadtxt('timesteps.txt',dtype='str')
i = 21 # final step in sim

#for i in range(len(timesteps)):
GM4_sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00'+timesteps[i])

h = GM4_sim.halos()
pynbody.analysis.halo.center(h[1],mode='ssc')
h1 = h[1]

pynbody.analysis.angmom.faceon(h1)
GM4_sim.physical_units()

# We want to track particles that NEVER become stars in GM4
madestarsinGM4_mask = np.in1d(GM1h1z0igasorder,GM4_sim.s['igasorder'])
GM1h1z0igasorder = GM1h1z0igasorder[~madestarsinGM4_mask]   # matched gas that never made stars
print(len(GM1h1z0igasorder))
GM1h1z0igasorder = np.random.choice(GM1h1z0igasorder,size=2)
#print(GM1h1z0igasorder)

GM4_GM1h1z0_gasmask  = np.in1d(GM4_sim.g['iord'],GM1h1z0igasorder,assume_unique=False)
#GM4_GM1h1z0_starmask = np.in1d(GM4_sim.s['igasorder'],GM1h1z0iords,assume_unique=False) 
print(len(GM4_sim.g['iord'][GM4_GM1h1z0_gasmask]))

gas_iord = GM4_sim.g['iord'][GM4_GM1h1z0_gasmask]
gas_mass = GM4_sim.g['mass'][GM4_GM1h1z0_gasmask]
gas_pos  = np.array(GM4_sim.g['pos'])[GM4_GM1h1z0_gasmask]
gas_v    = np.array(GM4_sim.g['vel'])[GM4_GM1h1z0_gasmask]
gas_rho  = np.array(GM4_sim.g['rho'])[GM4_GM1h1z0_gasmask]
gas_T    = GM4_sim.g['temp'][GM4_GM1h1z0_gasmask]
gas_met  = GM4_sim.g['metals'][GM4_GM1h1z0_gasmask]

print(GM4_sim.g['mass'].units,GM4_sim.g['pos'].units,GM4_sim.g['vel'].units,GM4_sim.g['rho'].units,GM4_sim.g['temp'].units)

import csv
with open('GM4_ulysses_data.csv', 'w') as csvfile:
    fieldnames = ['iords', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'rho', 'temp', 'metals']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=',')



    writer.writeheader()
    writer.writerow({'iords': str("%.3f" % np.max(gas_iord)), 'mass': str("%.3f" % np.max(gas_mass)), 'x': str("%.3f" % np.max(gas_pos)), 'y': str("%.3f" % np.max(gas_pos)), 'z': str("%.3f" % np.max(gas_pos)), 'vx': str("%.3f" % np.max(gas_v)), 'vy': str("%.3f" % np.max(gas_v)), 'vz': str("%.3f" % np.max(gas_v)), 'rho': str("%.3f" % np.max(gas_rho)), 'temp': str("%.3f" % np.max(gas_T)), 'metals': str("%.3f" % np.max(gas_met))})
    writer.writerow({'iords': str("%.3f" % np.min(gas_iord)), 'mass': str("%.3f" % np.min(gas_mass)), 'x': str("%.3f" % np.min(gas_pos)), 'y': str("%.3f" % np.min(gas_pos)), 'z': str("%.3f" % np.min(gas_pos)), 'vx': str("%.3f" % np.min(gas_v)), 'vy': str("%.3f" % np.min(gas_v)), 'vz': str("%.3f" % np.min(gas_v)), 'rho': str("%.3f" % np.min(gas_rho)), 'temp': str("%.3f" % np.min(gas_T)), 'metals': str("%.3f" % np.min(gas_met))})

    for i in range(len(gas_iord)):
        #print(str(gas_iord[i]), str(gas_mass[i]), str(gas_pos[i,0]))
        writer.writerow({'iords': str("%.3f" % gas_iord[i]), 'mass': str("%.3f" % gas_mass[i]), 'x': str("%.3f" % gas_pos[i,0]), 'y': str("%.3f" % gas_pos[i,1]), 'z': str("%.3f" % gas_pos[i,2]), 'vx': str("%.3f" % gas_v[i,0]), 'vy': str("%.3f" % gas_v[i,1]), 'vz': str("%.3f" % gas_v[i,2]), 'rho': str("%.3f" % gas_rho[i]), 'temp': str("%.3f" % gas_T[i]), 'metals': str("%.3f" % gas_met[i])})


np.savetxt('ulysses_gasiords.txt',gas_iord)    
