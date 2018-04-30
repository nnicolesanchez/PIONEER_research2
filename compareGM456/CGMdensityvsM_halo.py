# This script creates a plot of the CGM density 
# as a function of the MAIN HALO mass in the PIONEER
# suite of ChaNGa Nbody Simulations:
# P4, GM5, GM6 
# Middle main halo mass, more, and less
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

plt.rc('font', size=12, family='serif', style='normal', variant='normal', stretch='normal', weight='normal')
plt.rc('xtick', labelsize=12)
plt.rc('xtick.major', size=6, width=1)
plt.rc('lines', lw=2)
plt.rc('axes', lw=1, labelsize=12)

#P0 = pynbody.load('/nobackupp8/fgoverna/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.004096')             # Satellite Halo ID: h9
#GM1 = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.004096')     # Satellite Halo ID: h10
#GM4 = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.004096')    # Satellite Halo ID: h4
#GM5 = pynbody.load('/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.004096')
#GM6 = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.004096')


sims = ['/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.004096','/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.004096','/nobackupp8/fgoverna/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.004096']
names  = ['GM4','GM5','GM6']
colors = ['DodgerBlue','SteelBlue','IndianRed']

#timesteps_P0 = ['0136','0186','0223','0273','0345','0384','0454','0635','0768','0972','1152','1536','1739','1920','2304','2554','2688','3072','3195','3456','3840','4096']
#timesteps_GM1 = [0128,0136,0186,0192,0256,0273,0320,0384,0448,0454,0512,0576,0640,0704,0768,0832,0896,0960,0972,1024,1088,1152,1216,1280,1344,1408,1472,1536,1600,1664,1728,1739,1792,1856,1920,1984,2048,2112,2176,2240,2304,2368,2432,2496,2554,2560,2624,2688,2752,2816,2880,2944,3008,3072,3136,3200,3264,3328,3392,3456,3520,3584,3648,3712,3776,3840,3904,3968,4032,4096]
#timesteps_GM4 = [0128,0136,0186,0223,0256,0273,0345,0384,0454,0512,0635,0640,0768,0896,0972,1024,1152,1280,1408,1536,1664,1739,1792,1920,2048,2176,2304,2432,2554,2560,2688,2816,2944,3072,3195,3200,3328,3456,3584,3712,3840,3968,4096]

for i in range(len(sims)):
    sim = pynbody.load(sims[i])

    print('simulations at z = ','%.2f' % sim.properties['z'] )

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
    
    cool_stuff = (CGM_gas['temp'] < 10**5) & (CGM_gas['temp'] > 10**4)
    cool_gas = CGM_gas[cool_stuff]
    
    total_cool_mass = np.sum(cool_gas['mass'])
    print('Total mass of gas between 10^4-10^5 K:',total_cool_mass)
    
    
    total_halo_mass = h1['mass'].sum().in_units('Msol')
    print('Total mass of main halo '+str(h1['mass'].sum().in_units('Msol'))+' via pynbody')
    
    plt.plot(total_cool_mass,total_halo_mass,label=names[i],marker='*')

plt.xlabel(r'Total CGM Mass [M$_{sol}$] between 10$^4$-10$^5$ K')
plt.ylabel('Main Halo Mass [M$_{sol}$]')
plt.legend()
plt.savefig('H1massvsCGMcoolgasmass.pdf')
plt.show()
