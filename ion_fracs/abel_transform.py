import matplotlib.pyplot as plt
import numpy as np
from hdf5_Novi_R import *
import pynbody

m_p = 1.6726 * 10**-24

P0 = '/nobackupp2/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.003456'
GM1 = '/nobackupp2/nnsanche/pioneer50h243GM1.1536gst1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.003456'
GM2 = '/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.003456'
GM3 = '/nobackupp2/nnsanche/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.003456'

P0s = pynbody.load(P0)
GM1s = pynbody.load(GM1)
GM2s = pynbody.load(GM2)
GM3s = pynbody.load(GM3)

P0h = P0s.halos()
P0_h1 = P0h[1]
pynbody.analysis.angmom.faceon(P0_h1)
P0_CGMgas = P0_h1.g[P0_h1.g['r'].in_units('kpc') > 10]
P0_CGMgas_oden = P0_CGMgas['rho'].in_units('g cm^-3')*P0_CGMgas['OxMassFrac']/(16*1.6726 * 10**-24)
P0_CGMgas['ovi'] = hdf5_ion_frac(P0_CGMgas,ion='ovi')

GM1h = GM1s.halos()
GM1_h1 = GM1h[1]
pynbody.analysis.angmom.faceon(GM1_h1)
GM1_CGMgas = GM1_h1.g[GM1_h1.g['r'].in_units('kpc') > 10]
GM1_CGMgas_oden = GM1_CGMgas['rho'].in_units('g cm^-3')*GM1_CGMgas['OxMassFrac']/(16*1.6726 * 10**-24)
GM1_CGMgas['ovi'] = hdf5_ion_frac(GM1_CGMgas,ion='ovi')

GM2h = GM2s.halos()
GM2_h1 = GM2h[1]
pynbody.analysis.angmom.faceon(GM2_h1)
GM2_CGMgas = GM2_h1.g[GM2_h1.g['r'].in_units('kpc') > 10]
GM2_CGMgas_oden = GM2_CGMgas['rho'].in_units('g cm^-3')*GM2_CGMgas['OxMassFrac']/(16*1.6726 * 10**-24)
GM2_CGMgas['ovi'] = hdf5_ion_frac(GM2_CGMgas,ion='ovi')

GM3h = GM3s.halos()
GM3_h1 = GM3h[1]
pynbody.analysis.angmom.faceon(GM3_h1)
GM3_CGMgas = GM3_h1.g[GM3_h1.g['r'].in_units('kpc') > 10]
GM3_CGMgas_oden = GM3_CGMgas['rho'].in_units('g cm^-3')*GM3_CGMgas['OxMassFrac']/(16*1.6726 * 10**-24)
GM3_CGMgas['ovi'] = hdf5_ion_frac(GM3_CGMgas,ion='ovi')

b_kpc = np.linspace(10,260,num=26)
b_kpc_halfstep = np.linspace(5,255,num=26)
GM3_nOVI_kpc = np.interp(b_kpc,GM3_CGMgas['r'].in_units('kpc'),GM3_CGMgas['ovi']*GM3_CGMgas_oden)
GM3_N_OVI = np.zeros(len(b_kpc))
for i in range(len(b_kpc)):
    for j in range(len(b_kpc)):
        if((b_kpc[j]>=b_kpc_halfstep[i]) & (b_kpc[j]<b_kpc_halfstep[-1])):
            GM3_path = (b_kpc_halfstep[j+1]-b_kpc_halfstep[j])*b_kpc_halfstep[j+1]/np.sqrt(b_kpc_halfstep[j+1]**2-b_kpc[i]**2)
            GM3_N_OVI[i] += 2.*GM3_path*GM3_nOVI_kpc[j]*3.086*10**21
print(GM3_N_OVI)
import matplotlib.pyplot as plt
plt.plot(b_kpc,GM3_N_OVI)
b_kpc = np.linspace(10,260,num=26)
b_kpc_halfstep = np.linspace(5,255,num=26)
GM2_nOVI_kpc = np.interp(b_kpc,GM2_CGMgas['r'].in_units('kpc'),GM2_CGMgas['ovi']*GM2_CGMgas_oden)
GM2_N_OVI = np.zeros(len(b_kpc))
for i in range(len(b_kpc)):
    for j in range(len(b_kpc)):
        if((b_kpc[j]>=b_kpc_halfstep[i]) & (b_kpc[j]<b_kpc_halfstep[-1])):
            GM2_path = (b_kpc_halfstep[j+1]-b_kpc_halfstep[j])*b_kpc_halfstep[j+1]/np.sqrt(b_kpc_halfstep[j+1]**2-b_kpc[i]**2)
            GM2_N_OVI[i] += 2.*GM2_path*GM2_nOVI_kpc[j]*3.086*10**21
plt.plot(b_kpc,GM2_N_OVI)
b_kpc = np.linspace(10,260,num=26)
b_kpc_halfstep = np.linspace(5,255,num=26)
GM1_nOVI_kpc = np.interp(b_kpc,GM1_CGMgas['r'].in_units('kpc'),GM1_CGMgas['ovi']*GM1_CGMgas_oden)
GM1_N_OVI = np.zeros(len(b_kpc))
for i in range(len(b_kpc)):
    for j in range(len(b_kpc)):
        if((b_kpc[j]>=b_kpc_halfstep[i]) & (b_kpc[j]<b_kpc_halfstep[-1])):
            GM1_path = (b_kpc_halfstep[j+1]-b_kpc_halfstep[j])*b_kpc_halfstep[j+1]/np.sqrt(b_kpc_halfstep[j+1]**2-b_kpc[i]**2)
            GM1_N_OVI[i] += 2.*GM1_path*GM1_nOVI_kpc[j]*3.086*10**21
plt.plot(b_kpc,GM1_N_OVI)
b_kpc = np.linspace(10,260,num=26)
b_kpc_halfstep = np.linspace(5,255,num=26)
P0_nOVI_kpc = np.interp(b_kpc,P0_CGMgas['r'].in_units('kpc'),P0_CGMgas['ovi']*P0_CGMgas_oden)
P0_N_OVI = np.zeros(len(b_kpc))
for i in range(len(b_kpc)):
    for j in range(len(b_kpc)):
        if((b_kpc[j]>=b_kpc_halfstep[i]) & (b_kpc[j]<b_kpc_halfstep[-1])):
            P0_path = (b_kpc_halfstep[j+1]-b_kpc_halfstep[j])*b_kpc_halfstep[j+1]/np.sqrt(b_kpc_halfstep[j+1]**2-b_kpc[i]**2)
            P0_N_OVI[i] += 2.*GM1_path*GM1_nOVI_kpc[j]*3.086*10**21

P0_CGMprofile = profile.Profile(P0_CGMgas,min='0.1 kpc',max='250 kpc')
GM1_CGMprofile = profile.Profile(GM1_CGMgas,min='0.1 kpc',max='250 kpc')
GM2_CGMprofile = profile.Profile(GM2_CGMgas,min='0.1 kpc',max='250 kpc')
GM3_CGMprofile = profile.Profile(GM3_CGMgas,min='0.1 kpc',max='250 kpc')

plt.plot(GM2_CGMprofile['rbins'].in_units('kpc'),np.log10(GM2_CGMprofile['mass'].in_units('g')*GM2_CGMprofile['OxMassFrac']*GM2_CGMprofile['ovi']/(16*m_p)/GM2_CGMprofile._binsize.in_units('cm**2')),color='SteelBlue',linestyle='-')
plt.plot(GM3_CGMprofile['rbins'].in_units('kpc'),np.log10(GM3_CGMprofile['mass'].in_units('g')*GM3_CGMprofile['OxMassFrac']*GM3_CGMprofile['ovi']/(16*m_p)/GM3_CGMprofile._binsize.in_units('cm**2')),color='SkyBlue',linestyle='-')
plt.plot(GM1_CGMprofile['rbins'].in_units('kpc'),np.log10(GM1_CGMprofile['mass'].in_units('g')*GM1_CGMprofile['OxMassFrac']*GM1_CGMprofile['ovi']/(16*m_p)/GM1_CGMprofile._binsize.in_units('cm**2')),color='FireBrick',linestyle='-')
plt.plot(P0_CGMprofile['rbins'].in_units('kpc'),np.log10(P0_CGMprofile['mass'].in_units('g')*P0_CGMprofile['OxMassFrac']*P0_CGMprofile['ovi']/(16*m_p)/P0_CGMprofile._binsize.in_units('cm**2')),color='Red',linestyle='-')
plt.plot(b_kpc,np.log10(P0_N_OVI),color='SteelBlue',linestyle='--')
plt.plot(b_kpc,np.log10(GM1_N_OVI),color='SkyBlue',linestyle='--')
plt.plot(b_kpc,np.log10(GM2_N_OVI),color='FireBrick',linestyle='--')
plt.plot(b_kpc,np.log10(GM3_N_OVI),color='Red',linestyle='--')
plt.ylabel(r'N$_{OVI}$')
plt.xlabel('R/kpc')
plt.show()
