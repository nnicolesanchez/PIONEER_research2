# This script creates a plot of the fraction of OVI
# gas in the simulations GM1 & GM3(GM4 in old terms)
# at a variety of densities

# N. Nicole Sanchez -- Created: July 9, 2018
# Univ. of Wash, S  -- Edited:
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from hdf5_Novi_R import *
import numpy as np
import pynbody
import sys

if len(sys.argv) == 1:
    print('No galaxy selected. Current options: P0, GM1, GM7, GM4')
    print('Syntax: "f_ovi_vs_T.py GM1"')
    quit()
else:
    if (str(sys.argv[1]) == 'P0'):
        sim = pynbody.load('/nobackupp2/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.003456')
    elif (str(sys.argv[1]) == 'GM1'):
        sim = pynbody.load('/nobackupp2/nnsanche/pioneer50h243GM1.1536gst1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.003456')
    elif (str(sys.argv[1]) == 'GM4'):
        sim = pynbody.load('/nobackupp2/nnsanche/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.003456')    
    elif (str(sys.argv[1]) == 'GM5'):
        sim = pynbody.load('/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.003546')
    elif (str(sys.argv[1]) == 'GM6'):
        sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.003456')
    elif (str(sys.argv[1]) == 'GM7'):
        sim = pynbody.load('/nobackupp2/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.003456') 

    else :
        print('Not a valid option. Current options: P0, GM1, GM7, GM4')
        print('Syntax: "f_ovi_vs_T.py GM1"')
        quit()    

    if str(sys.argv[1]) == 'GM4':
        name = 'GM3'
    elif str(sys.argv[1]) == 'GM7':
        name = 'GM2'
    else:
        name = str(sys.argv[1])
    print(name+' simulation at z = ','%.2f' % sim.properties['z'] )

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
z_max = 10   # kpc

Rg_d = ((h1.g['x'].in_units('kpc'))**2. + (h1.g['y'].in_units('kpc'))**2.)**(0.5)
disk_gas_xymax =  (Rg_d < r_max)
disk_gas_zmax  = (h1.g['z'].in_units('kpc') < z_max) & (h1.g['z'].in_units('kpc') > -z_max)

disk_gas_mask = disk_gas_xymax & disk_gas_zmax
disk_gas = h1.g[disk_gas_xymax & disk_gas_zmax]
CGM_gas  = h1.g[~disk_gas_mask]

density_array = [r'10$^{-2}$',r'10$^{-3}$',r'10$^{-4}$',r'10$^{-5}$']

# Using original .npz table from Stinson 2012
ovi_npz    = pynbody.analysis.ionfrac.calculate(CGM_gas,ion='ovi',mode='new') 

#rho_10_m2  = (CGM_gas['rho'] > 0.05) & (CGM_gas['rho'] > 0.005)
#CGMg_10_m2 = CGM_gas[rho_10_m2]

CGMg_fovi = (CGM_gas['mass']*CGM_gas['OxMassFrac']*ovi_npz)/CGM_gas['mass']
plt.hexbin(np.log10(CGM_gas['temp']),CGMg_fovi,C=np.log10(CGM_gas['rho']),cmap=cm.plasma)
#plt.plot(CGMg_10_m2['temp'],CGMg_10_m2['mass']*CGMg_10_m2['OxMassFrac']*ovi_npz[rho_10_m2],label=density_array[0],linestyle=None,marker='.')
plt.ylabel(r'f$_{OVI}$')
plt.xlabel('T [K]')
plt.colorbar(label=r'$\rho$')
plt.title(name+', using .npz table')
plt.savefig(name+'_fovi_T_npz.pdf')
plt.show()

# Using new HDF5 table from Corlies and Hummels
ovi_hdf5    = hdf5_ion_frac(CGM_gas,ion='ovi')

CGMg_fovi = (CGM_gas['mass']*CGM_gas['OxMassFrac']*ovi_hdf5)/CGM_gas['mass']
plt.hexbin(np.log10(CGM_gas['temp']),CGMg_fovi,C=np.log10(CGM_gas['rho']),cmap=cm.plasma)
plt.ylabel(r'f$_{OVI}$')
plt.xlabel('T [K]')
plt.colorbar(label=r'$\rho$')
plt.title(name+', using hdf5 table')
plt.savefig(name+'_fovi_T_hdf5.pdf')
plt.show()










