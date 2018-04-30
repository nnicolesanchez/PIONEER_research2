# This script creates a plot of the CGM density 
# as a function of the SATELLITE mass in the PIONEER
# suite of ChaNGa Nbody Simulations:
# P0, GM1, GM4 
# Most massive satellite mass, less, and less
# (since the thing that is changing bw these is SAT mass)

#     - Outputs:
#         1. plot of total mass of gas between 10^4 & 10^5 K 
#          as a function of satellite 
#         2. plot of CGM density as func of satellite


# N. Nicole Sanchez -- August 2017
# Univ. of Wash.    -- Nbody Shop
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np
import pynbody
import sys

plt.rc('font', size=12, family='serif', style='normal', variant='normal', stretch='normal', weight='no
rmal')
plt.rc('xtick', labelsize=12)
plt.rc('xtick.major', size=6, width=1)
plt.rc('lines', lw=2)
plt.rc('axes', lw=1, labelsize=12)

#if len(sys.argv) == 1:
#    print('No galaxy selected. Current options: P0, GM1, GM4, GM5, GM6')
#    print('Syntax: "GM_xygasplots.py GM1"')
#    quit()
#else:
#    if (str(sys.argv[1]) == 'P0'):
P0 = pynbody.load('/nobackupp8/fgoverna/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.004096')             # Satellite Halo ID: h9
#    elif (str(sys.argv[1]) == 'GM1'):
GM1 = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.004096')     # Satellite Halo ID: h10
#    elif (str(sys.argv[1]) == 'GM4'):
GM4 = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.004096')    # Satellite Halo ID: h4
#    elif (str(sys.argv[1]) == 'GM5'):
#        sim = pynbody.load('/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.004096')
#    elif (str(sys.argv[1]) == 'GM6'):
#        sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.004096')

#    else :
#        print('Not a valid option. Current options: P0, GM1, GM4')
#        print('Syntax: "GM_xygasplots.py GM1"')
#        quit()    

#    name = str(sys.argv[1])
print(name+' simulation at z = ','%.2f' % P0.properties['z'] )

sim.properties
sim.properties['time'].in_units('Gyr')
sim.loadable_keys()
sim.physical_units()
h = sim.halos()
h1 = h[1]
pynbody.analysis.angmom.faceon(h1)

# Constants
m_H = 1.6733 * 10**-24 #g

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

