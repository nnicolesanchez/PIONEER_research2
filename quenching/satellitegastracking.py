# This script reads in the satellite's particles before the merger 
# and sees where they go:
#     - Do they settle in the disk? 
#     - How do I include gas that has come in with the merging galaxy?

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import os.path
import seaborn as sns
import pandas as pd
import numpy as np
import pynbody
from pynbody.analysis import profile

sim = ['/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.00']
labels = ['P0','GM1','GM4','GM5','GM6','GM7']
pre_merger_time = '1739'
z_0             = '4096' 

























@pynbody.derived_array
def N_OVI(f):
    ovi = pynbody.analysis.ionfrac.calculate(f,ion='ovi',mode='new')
    m_p = 1.6726 * 10**-24 #g
    return f.gas['rho'].in_units('g cm**-3')*ovi*f.gas['OxMassFrac']/(16*m_p)



ts = np.loadtxt('../'+labels[k]+'/timesteps.txt',dtype=str)
t = len(ts)-1
for t in range(6,len(ts)):
    print('Loading sim:',sim[k],' at timestep:',ts[t])
    if t == len(ts)-1 :
        continue

    ######################
    # READ IN SIMULATION #
    ######################
    f = pynbody.load(sim[k]+ts[t])
    f.physical_units()
    h = f.halos()
    h1 = h[1]
    pynbody.analysis.halo.center(h1,mode='ssc')
    pynbody.analysis.angmom.faceon(h1)

    ###################
    # ISOLATE CGM GAS #   
    ###################
    # Isolate and remove disk stars within radius 0-10 kpc & vertically 10 kpc 
    r_max = 10  # kpc
    twenty_kpc_incm = 6.171*(10**22)
    
    Rg_d = ((h1.g['x'].in_units('kpc'))**2. + (h1.g['y'].in_units('kpc'))**2. + (h1.g['z'].in_units('kpc'))**2.)**(0.5)
    print('Max radial gas distance',np.max(Rg_d))

    if (np.max(Rg_d) < 50):
        print('Halo is too small for gas at 50 kpc; skipping this step')
        continue

    disk_gas_xyzmax =  (Rg_d < r_max)
    disk_gas_mask = disk_gas_xyzmax #& disk_gas_zmax
    disk_gas = h1.g[disk_gas_mask] #& disk_gas_zmax]
    CGM_gas  = h1.g[~disk_gas_mask]
    CGM_temp = np.array(CGM_gas['temp'])
    
    #########################
    # CALCULATE OVI DENSITY #
    #########################
    # Ionization fraction of OVI compare to total Oxygen
    CGM_gas['ovi'] = pynbody.analysis.ionfrac.calculate(CGM_gas,ion='ovi',mode='new') 
    m_p = 1.6726 * 10**-24 #g
    
    print('Total mass in CGM:', np.sum(CGM_gas['mass']))
    #print('Total mass in metals:',np.sum(CGM_gas['mass']*CGM_gas['metals'])) # 'metals' *IS* metallicity
    #print('Total mass Oxygen in CGM:', np.sum(CGM_gas['OxMassFrac']*CGM_gas['mass']),CGM_gas['mass'].units)
    #print('Total mass in OVI in CGM:', np.sum(CGM_gas['OxMassFrac']*CGM_gas['mass']*CGM_gas['ovi']))
    #print('OVI fractions:',np.average(CGM_gas['ovi']))

    from pynbody.analysis import profile
    p = profile.Profile(CGM_gas,min='0.1 kpc',max='250 kpc')

#    print(p['rbins'].in_units('kpc'),(p['mass'].in_units('g')*p['OxMassFrac']*p['ovi']/(16*m_p))/p._binsize.in_units('cm**2'))
    mask_50kpc = (np.abs(p['rbins'].in_units('kpc')-50)).argmin()
#    print('50 kpc bin',p['rbins'].in_units('kpc')[mask_50kpc],(p['mass'].in_units('g')[mask_50kpc]*p['OxMassFrac']*p['ovi']/(16*m_p))/p._binsize.in_units('cm**2'))

    rbin_50kpc_thru_time.append(p['rbins'].in_units('kpc')[mask_50kpc])
    N_ovi_50kpc_thru_time.append((p['mass'].in_units('g')[mask_50kpc]*p['OxMassFrac'][mask_50kpc]*p['ovi'][mask_50kpc]/(16*m_p))/p._binsize.in_units('cm**2')[mask_50kpc])
    print(N_ovi_50kpc_thru_time)
    times.append(f.properties['time'].in_units('Gyr'))




