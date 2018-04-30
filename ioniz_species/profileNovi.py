import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import os.path
import seaborn as sns
import pandas as pd
import numpy as np
import pynbody
from pynbody.analysis import profile

@pynbody.derived_array
def N_OVI(f):
    ovi = pynbody.analysis.ionfrac.calculate(f,ion='ovi',mode='new')
    m_p = 1.6726 * 10**-24 #g
    return f.gas['rho'].in_units('g cm**-3')*ovi*f.gas['OxMassFrac']/(16*m_p)

# Just using k = 1 and k = 2, for GM1 & GM4 for now
k = 0
## MOVED FILES FROM FABIO TO ALYSON BROOKS: /nobackupp8/ambrook2/fgoverna_pleiades_p8_files
sim = ['/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.00']
labels = ['P0','GM1','GM4','GM5','GM6','GM7']
colors = sns.cubehelix_palette(8)
print('LOADING SIM:',labels[k])

ts = np.loadtxt('../'+labels[k]+'/timesteps.txt',dtype=str)
t = len(ts)-1
if k == 5:
    t = len(ts)-2
#for t in range(len(ts)):
print('Loading sim:',sim[k],' at timestep:',ts[t])

######################
# READ IN SIMULATION #
######################
f = pynbody.load(sim[k]+ts[t])
pynbody.analysis.halo.center(f.star)
f.physical_units()
h = f.halos()
h1 = h[1]
pynbody.analysis.angmom.faceon(h1)

###################
# ISOLATE CGM GAS #   
###################
# Isolate and remove disk stars within radius 0-10 kpc & vertically 10 kpc 
r_max = 10  # kpc
twenty_kpc_incm = 6.171*(10**22)

Rg_d = ((h1.g['x'].in_units('kpc'))**2. + (h1.g['y'].in_units('kpc'))**2. + (h1.g['z'].in_units('kpc'))**2.)**(0.5)
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


COS = pd.read_csv('COShalo_obs.txt',header=0,delim_whitespace=True,comment='#',index_col=False)
COS_ID = COS['ID']
COS_OVI = COS['LogN0VI']
COS_Rkpc = COS['Rho(kpc)'][COS_OVI != 0]
COS_OVI = COS_OVI[COS_OVI != 0]
#        print(COS)
plt.plot(COS_Rkpc[COS['RED?'] == 'yes'],COS_OVI[COS['RED?'] == 'yes'],marker='.',color='Red',linestyle=' ',label='COS Ellipticals')
plt.plot(COS_Rkpc[COS['RED?'] == 'no'],COS_OVI[COS['RED?'] == 'no'],marker='.',color='DodgerBlue',linestyle=' ',label='COS Spirals')

profile = profile.Profile(CGM_gas,min='0.1 kpc',max='250 kpc')

#print('Bins',profile['bins'])
#print('RBins',profile['dr'])
#print(profile['mass'].units)
#plt.plot(profile['rbins'].in_units('kpc'),(profile['mass'].in_units('g cm**-3')*profile['OxMassFrac']*profile['ovi']/(16*m_p))/profile.bins,'r-')
#plt.plot(profile['rbins'].in_units('kpc'),(profile['rho'].in_units('g cm**-3')*profile['rbins'].in_units('cm'),'r-'))


np.savetxt(labels[k]+'_N_OVI_z0.np',np.log10((profile['mass'].in_units('g')*profile['OxMassFrac']*profile['ovi']/(16*m_p))/profile._binsize.in_units('cm**2')))
np.savetxt(labels[k]+'_T_z0.np',np.log10(profile['temp'].in_units('K')))
np.savetxt(labels[k]+'_rho_z0.np',np.log10(profile['rho'].in_units('g cm^-3')))
np.savetxt(labels[k]+'_Rbins_z0.np',profile['rbins'].in_units('kpc'))
quit()
plt.plot(profile['rbins'].in_units('kpc'),np.log10((profile['mass'].in_units('g')*profile['OxMassFrac']*profile['ovi']/(16*m_p))/profile._binsize.in_units('cm**2')),label=labels[k])
plt.title('z = 0.00')
plt.ylabel(r'N$_{OVI}$ (cm$^{-2}$)')
plt.xlabel('R (kpc)')
plt.ylim(13,16)
plt.xlim(-10,260)
plt.legend()
plt.savefig('Novi_profile_'+labels[k]+'_AHF.pdf')
plt.show()


plt.plot(profile['rbins'].in_units('kpc'),np.log10(profile['temp'].in_units('K')))
plt.ylabel('Temp [K]')
plt.xlabel('R (kpc)')
plt.savefig('T_profile_'+labels[k]+'.pdf')
plt.show()
