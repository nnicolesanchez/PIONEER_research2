# This script creates the column density plots for the GM
# suite for h243 (P0, GM1, GM2(GM7), GM3(GM4))
# With and without black hole physics

# This script utilizes the hdf5 table of ion fractions 
# from Trident (from Corlies & Hummels)

# N. Nicole Sanchez  -- Created: June 26, 2018
# Univ of W, Seattle -- Edited: June 26, 2018
from pynbody.analysis.interpolate import _interpolate3d
from pynbody.analysis import profile
import matplotlib.pyplot as plt
from hdf5_Novi_R import *
import pandas as pd
import numpy as np
import pynbody
import h5py
import sys

if len(sys.argv) == 1:
    print('No galaxy selected. Current options: P0, GM1, GM7, GM4')
    print('For galaxies without BH physics: "P0_noBH"')
    print('Syntax: "make_profiles_3456.py GM1"')
    quit()
elif (str(sys.argv[1]) == 'P0'):
    lab = str(sys.argv[1])
    i = 0
elif (str(sys.argv[1]) == 'GM1'):
    lab = str(sys.argv[1])
    i = 1
elif (str(sys.argv[1]) == 'GM7'):
    lab = str(sys.argv[1])
    i = 2
elif (str(sys.argv[1]) == 'GM4'):
    lab = str(sys.argv[1])
    i = 3
elif (str(sys.argv[1]) == 'P0noBH'):
    lab = str(sys.argv[1])
    i = 4
elif (str(sys.argv[1]) == 'GM1noBH'):
    lab = str(sys.argv[1])
    i = 5
elif (str(sys.argv[1]) == 'GM7noBH'):
    lab = str(sys.argv[1])
    i = 6
elif (str(sys.argv[1]) == 'GM4noBH'):
    lab = str(sys.argv[1])
    i = 7
else:
    print('Invalid option. Goodbye.')
    quit()

sims      = ['/nobackupp2/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.003456','/nobackupp2/nnsanche/pioneer50h243GM1.1536gst1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.003456','/nobackupp2/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.003456','/nobackupp2/nnsanche/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.003456','/nobackup/nnsanche/NO_BHs/pioneer50h243.1536gst1bwK1_3456/pioneer50h243.1536gst1bwK1.003456','/nobackup/nnsanche/NO_BHs/pioneer50h243GM1.1536gst1bwK1_3456/pioneer50h243GM1.1536gst1bwK1.003456','/nobackup/nnsanche/NO_BHs/pioneer50h243GM7.1536gst1bwK1_3456/pioneer50h243GM7.1536gst1bwK1.003456','/nobackup/nnsanche/NO_BHs/pioneer50h243GM4.1536gst1bwK1_3456/pioneer50h243GM4.1536gst1bwK1.003456']

ion_labels = ['oi','oii','oiii','oiv','ov','ovi','ovii','oviii']
nion       = [1,2,3,4,5,6,7,8]
m_p = 1.6726 * 10**-24 #g

print(lab)
f = pynbody.load(sims[i])
pynbody.analysis.halo.center(f.star)
f.physical_units()
h = f.halos()
h1 = h[1]
#pynbody.analysis.angmom.faceon(h1)
pynbody.analysis.angmom.sideon(h1)

#r_max = 10  # kpc
#twenty_kpc_incm = 6.171*(10**22)
#Rg_d = ((h1.g['x'].in_units('kpc'))**2. + (h1.g['y'].in_units('kpc'))**2. + (h1.g['z'].in_units('kpc'))**2.)**(0.5)
#disk_gas_xyzmax =  (Rg_d < r_max)
#disk_gas_mask = disk_gas_xyzmax #& disk_gas_zmax
#disk_gas = h1.g[disk_gas_mask] #& disk_gas_zmax]

CGM_gas  = h1.g[h1.g['r'].in_units('kpc') < 10]
CGM_temp = np.array(CGM_gas['temp'])

for j in range(len(ion_labels)):
    CGM_gas[ion_labels[j]] = hdf5_ion_frac(CGM_gas,ion=ion_labels[j]) 
   
CGMprofile = profile.Profile(CGM_gas,min='10 kpc',max='250 kpc')

for k in range(len(ion_labels)):
    np.savetxt(lab+'/'+lab+'_N'+ion_labels[k]+'_3456.np',np.log10((CGMprofile['mass'].in_units('g')*CGMprofile['OxMassFrac']*CGMprofile[ion_labels[k]]/(16*m_p))/CGMprofile._binsize.in_units('cm**2')))
    np.savetxt(lab+'/'+lab+'_'+ion_labels[k]+'frac_3456.np',CGMprofile[ion_labels[k]])

np.savetxt(lab+'/'+lab+'_T_3456.np',np.log10(CGMprofile['temp'].in_units('K')))
np.savetxt(lab+'/'+lab+'_rho_3456.np',np.log10(CGMprofile['rho'].in_units('g cm^-3')))
np.savetxt(lab+'/'+lab+'_Rbins_3456.np',CGMprofile['rbins'].in_units('kpc'))
np.savetxt(lab+'/'+lab+'_totgasmass_3456.np',CGMprofile['mass'].in_units('Msol'))
np.savetxt(lab+'/'+lab+'_Z_3456.np',CGMprofile['metals'])
np.savetxt(lab+'/'+lab+'_Omass_3456.np',CGMprofile['OxMassFrac']*CGMprofile['mass'].in_units('Msol')) # in units Msol
