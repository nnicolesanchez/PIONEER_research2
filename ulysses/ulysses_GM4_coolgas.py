import matplotlib.pyplot as plt
#import pandas as pd
import numpy as np
import pynbody

N = 2000

#sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.004096')
sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.004096')

sim.properties
sim.properties['time'].in_units('Gyr')
sim.loadable_keys()
sim.physical_units()
h = sim.halos()
h1 = h[1]
pynbody.analysis.angmom.faceon(h1)
pynbody.analysis.halo.center(h[1],mode='ssc')

print('simulations at z = ','%.2f' % sim.properties['z'] )

# Want to isolate CGM  
# Isolate and remove disk stars within radius 0-10 kpc & vertically 4 kpc 
r_max = 10  # kpc
z_max = 4   # kpc

Rg_d = ((h1.g['x'].in_units('kpc'))**2. + (h1.g['y'].in_units('kpc'))**2.)**(0.5)
disk_gas_xymax =  (Rg_d < r_max)
disk_gas_zmax  = (h1.g['z'].in_units('kpc') < z_max) & (h1.g['z'].in_units('kpc') > -z_max)
    
disk_gas_mask = disk_gas_xymax & disk_gas_zmax
disk_gas = h1.g[disk_gas_xymax & disk_gas_zmax]
CGM_gas  = h1.g[~disk_gas_mask]
    
cool_stuff = (CGM_gas['temp'] < 10**5) & (CGM_gas['temp'] > 10**4.5)
cool_gas = CGM_gas[cool_stuff]
print('Number of cool gas particles:',len(cool_gas))

pynbody.plot.sph.image(cool_gas,qty="rho",units="g cm^-3",width=100,cmap="Greys")
plt.savefig('GM1_coolgas_rho_faceon.pdf')
#plt.show()

iord_mask = np.random.choice(cool_gas['iord'],size=N,replace=False)

# TURN ON FOR IORDS PRIOR TO 4096
#iord_4096 = np.loadtxt('GM1_coolgas_ulysses_data_4096.csv',delimiter=',',skiprows=4,usecols=0)
#print(iord_4096)

mask = np.in1d(cool_gas['iord'],iord_mask,assume_unique=True)
print(len(cool_gas['iord'][mask]))

gas_iord = np.array(cool_gas['iord'][mask])
print('Number of cool gas particles:',len(cool_gas['iord'][mask]))
gas_mass = np.array(cool_gas['mass'][mask])
gas_pos  = np.array(cool_gas['pos'][mask])*sim.properties['a']
gas_v    = np.array(cool_gas['vel'][mask])*sim.properties['a']
gas_rho  = np.array(cool_gas['rho'][mask])
gas_T    = np.array(cool_gas['temp'][mask])
gas_met  = np.array(cool_gas['metals'][mask])

#pynbody.plot.sph.image(cool_gas[mask],qty="rho",units="g cm^-3",width=100,cmap="Greys")
#plt.show()


import csv
#with open('GM1_coolgas_ulysses_data_4096.csv', 'w') as csvfile:
with open('GM4_coolgas_ulysses_data_4096.csv', 'w') as csvfile:
    fieldnames = ['iords', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'rho', 'temp', 'metals']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=',')

    # Header
    writer.writeheader()

    # Units
#    writer.writerow({'iords': 'NoUnit()', 'mass': str(sim.g['mass'].units), 'x': str(sim.g['pos'].units)+'*a', 'y': str(sim.g['pos'].units)+'*a', 'z': str(sim.g['pos'].units)+'*a', 'vx': str(sim.g['vel'].units)+'*a', 'vy': str(sim.g['vel'].units)+'*a', 'vz': str(sim.g['vel'].units)+'*a', 'rho': 'log('+str(sim.g['rho'].units)+')','temp': 'log(NoUnit())', 'metals': 'log(NoUnit())'})

    # Max
#    writer.writerow({'iords': str(int(np.max(gas_iord))), 'mass': str("%.2f" % np.max(gas_mass)), 'x': str("%.2f" % np.max(gas_pos)), 'y': str("%.2f" % np.max(gas_pos)), 'z': str("%.2f" % np.max(gas_pos)), 'vx': str("%.2f" % np.max(gas_v)), 'vy': str("%.2f" % np.max(gas_v)), 'vz': str("%.2f" % np.max(gas_v)), 'rho': str("%.2f" % np.max(gas_rho)), 'temp': str("%.2f" % np.max(gas_T)), 'metals': str("%.2f" % np.max(gas_met))})
    
    # Min
#    writer.writerow({'iords': str(int(np.min(gas_iord))), 'mass': str("%.2f" % np.min(gas_mass)), 'x': str("%.2f" % np.min(gas_pos)), 'y': str("%.2f" % np.min(gas_pos)), 'z': str("%.2f" % np.min(gas_pos)), 'vx': str("%.2f" % np.min(gas_v)), 'vy': str("%.2f" % np.min(gas_v)), 'vz': str("%.2f" % np.min(gas_v)), 'rho': str("%.2f" % np.min(gas_rho)), 'temp': str("%.2f" % np.min(gas_T)), 'metals': str("%.2f" % np.min(gas_met))})


    for i in range(len(gas_iord)):
        print(str(int(gas_iord[i])), str(gas_mass[i]), str(gas_pos[i,0]))
        writer.writerow({'iords': str(gas_iord[i]), 'mass': str("%.2f" % gas_mass[i]), 'x': str("%.2f" % gas_pos[i,0]), 'y': str("%.2f" % gas_pos[i,1]), 'z': str("%.2f" % gas_pos[i,2]), 'vx': str("%.2f" % gas_v[i,0]), 'vy': str("%.2f" % gas_v[i,1]), 'vz': str("%.2f" % gas_v[i,2]), 'rho': str("%.2f" % np.log10(gas_rho[i])), 'temp': str("%.2f" % np.log10(gas_T[i])), 'metals': str("%.2f" % np.log10(gas_met[i]))})


