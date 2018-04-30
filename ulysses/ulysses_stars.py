cimport matplotlib.pyplot as plt
#import pandas as pd
import numpy as np
import pynbody

N = 100000

GM1h1z0igasorder = np.loadtxt('../compare_GM1_GM4/GM1_h1_tracedigasorder.txt')
print(len(GM1h1z0igasorder))

timesteps = np.loadtxt('timesteps.txt',dtype='str')
i = 21 # final step in sim

#for i in range(len(timesteps)):
GM1_sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00'+timesteps[i])

h = GM1_sim.halos()
pynbody.analysis.halo.center(h[1],mode='ssc')
h1 = h[1]

pynbody.analysis.angmom.faceon(h1)
GM1_sim.physical_units()

# We want to track star particles in GM1 that NEVER become stars in GM4
madestarsinGM4_mask = np.in1d(GM1h1z0igasorder,GM1_sim.s['igasorder'])
GM1h1z0igasorder = GM1h1z0igasorder[madestarsinGM4_mask]   # matched gas that made stars
print(len(GM1h1z0igasorder))
GM1h1z0igasorder = np.random.choice(GM1h1z0igasorder,size=N)
#print(GM1h1z0igasorder)

GM1h1z0_gasmask  = np.in1d(GM1_sim.s['igasorder'],GM1h1z0igasorder,assume_unique=False)
print(len(GM1_sim.s['igasorder'][GM1h1z0_gasmask]))

gas_iord = GM1_sim.s['igasorder'][GM1h1z0_gasmask]
gas_mass = GM1_sim.s['mass'][GM1h1z0_gasmask]
gas_pos  = np.array(GM1_sim.s['pos'])[GM1h1z0_gasmask]*GM1_sim.properties['a']
gas_v    = np.array(GM1_sim.s['vel'])[GM1h1z0_gasmask]*GM1_sim.properties['a']
gas_rho  = 0 #np.array(GM1_sim.s['rho'])[GM1h1z0_gasmask]
gas_T    = 0 # GM1_sim.s['temp'][GM1h1z0_gasmask]
gas_met  = GM1_sim.s['metals'][GM1h1z0_gasmask]

print(len(gas_iord))

import csv
with open('GM1_ulysses_data.csv', 'w') as csvfile:
    fieldnames = ['iords', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'rho', 'temp', 'metals']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=',')

    # Header
    writer.writeheader()

    # Units
    writer.writerow({'iords': 'NoUnit()', 'mass': str(GM1_sim.s['mass'].units), 'x': str(GM1_sim.s['pos'].units)+'*a', 'y': str(GM1_sim.s['pos'].units)+'*a', 'z': str(GM1_sim.s['pos'].units)+'*a', 'vx': str(GM1_sim.s['vel'].units)+'*a', 'vy': str(GM1_sim.s['vel'].units)+'*a', 'vz': str(GM1_sim.s['vel'].units)+'*a', 'rho': 'log('+str(GM1_sim.s['rho'].units)+')','temp': 'log(NoUnit())', 'metals': 'log(NoUnit())'})

    # Max
    writer.writerow({'iords': str(int(np.max(gas_iord))), 'mass': str("%.2f" % np.max(gas_mass)), 'x': str("%.2f" % np.max(gas_pos)), 'y': str("%.2f" % np.max(gas_pos)), 'z': str("%.2f" % np.max(gas_pos)), 'vx': str("%.2f" % np.max(gas_v)), 'vy': str("%.2f" % np.max(gas_v)), 'vz': str("%.2f" % np.max(gas_v)), 'rho': str("%.2f" % np.max(gas_rho)), 'temp': str("%.2f" % np.max(gas_T)), 'metals': str("%.2f" % np.max(gas_met))})
    
    # Min
    writer.writerow({'iords': str(int(np.min(gas_iord))), 'mass': str("%.2f" % np.min(gas_mass)), 'x': str("%.2f" % np.min(gas_pos)), 'y': str("%.2f" % np.min(gas_pos)), 'z': str("%.2f" % np.min(gas_pos)), 'vx': str("%.2f" % np.min(gas_v)), 'vy': str("%.2f" % np.min(gas_v)), 'vz': str("%.2f" % np.min(gas_v)), 'rho': str("%.2f" % np.min(gas_rho)), 'temp': str("%.2f" % np.min(gas_T)), 'metals': str("%.2f" % np.min(gas_met))})


    for i in range(len(gas_iord)):
        print(str(int(gas_iord[i])), str(gas_mass[i]), str(gas_pos[i,0]))
        writer.writerow({'iords': str(gas_iord[i]), 'mass': str("%.2f" % gas_mass[i]), 'x': str("%.2f" % gas_pos[i,0]), 'y': str("%.2f" % gas_pos[i,1]), 'z': str("%.2f" % gas_pos[i,2]), 'vx': str("%.2f" % gas_v[i,0]), 'vy': str("%.2f" % gas_v[i,1]), 'vz': str("%.2f" % gas_v[i,2]), 'rho': str("%.2f" % np.log10(gas_rho)), 'temp': str("%.2f" % np.log10(gas_T)), 'metals': str("%.2f" % np.log10(gas_met[i]))})


