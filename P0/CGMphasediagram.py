# This script creates a phase diagram (temp vs density)
# for the CGM of a ChaNGa Nbody Simulation

#     - Outputs:
#          plots of phase diagram


# N. Nicole Sanchez -- July 2 2017
# Univ. of Wash.    -- Nbody Shop
import matplotlib.pyplot as plt
import numpy as np
import pynbody
import sys

plt.rc('font', size=12, family='serif', style='normal', variant='normal', stretch='normal', weight='normal')
plt.rc('xtick', labelsize=12)
plt.rc('xtick.major', size=6, width=1)
plt.rc('lines', lw=2)
plt.rc('axes', lw=1, labelsize=12)

#if len(sys.argv) == 1:
#    print('Syntax: python patient0_scalelength.py max_age min_age r_max z_max')
#    age_max = 6.
#else:
#    age_max = float(sys.argv[1])

timesteps = np.loadtxt('../Patient0/timesteps.txt',dtype='str')

i = 21
#for i in range(len(timesteps)):
print('Starting timestep: ', timesteps[i])
sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.002304')
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

Rg_d = ((h1.g['x'].in_units('kpc'))**2. + (h1.g['y'].in_units('kpc'))**2.)**(0.5)
disk_gas_xymax =  (Rg_d < r_max)
disk_gas_zmax  = (h1.g['z'].in_units('kpc') < z_max) & (h1.g['z'].in_units('kpc') > -z_max)

disk_gas_mask = disk_gas_xymax & disk_gas_zmax
disk_gas = h1.g[disk_gas_xymax & disk_gas_zmax]
CGM_gas  = h1.g[~disk_gas_mask]
mini_CGM_gas = np.random.choice(CGM_gas['iord'],size=10000)
mini_CGM_gas_mask = np.in1d(CGM_gas['iord'],mini_CGM_gas)

plt.scatter(np.log10(CGM_gas['rho'][mini_CGM_gas_mask]),np.log10(CGM_gas['temp'][mini_CGM_gas_mask]),c=np.log10(CGM_gas['metals'][mini_CGM_gas_mask]),cmap=plt.get_cmap('jet'), alpha=0.5,s=np.log10(CGM_gas['mass'][mini_CGM_gas_mask]))
plt.ylabel(r'log$_{10}$(T ('+str(CGM_gas['temp'].units)+'))')
plt.xlabel(r'log$_{10}$($\rho$ (M$_{\odot}$ kpc$^{-3}$))')
plt.colorbar(label=(r'log$_{10}$(metal fraction)'))
plt.title('P0 at z=0')
plt.savefig('P0_phasediagram.pdf')
plt.show()
