import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pynbody.analysis import profile
import pynbody

labels = ['P0','GM1','GM4','GM5','GM6','GM7']
colors = ['DodgerBlue','SteelBlue','FireBrick','IndianRed','Salmon','Orange']

k = 2
ts = np.loadtxt('../'+labels[k]+'/timesteps.txt',dtype=str)
z  = np.loadtxt('../'+labels[k]+'/redshifts.txt',dtype=float,delimiter=',',usecols=1)
ts = np.loadtxt('../'+labels[k]+'/timesteps.txt',dtype=str)

for j in range(len(labels)):
    #t = len(ts)-1
    #if j == 5:
    #    t = len(ts)-2
    t = len(ts)-6
    print('Timestep: ',ts[t],'Redshift:','%.2f' % z[t])

    Novi = np.loadtxt('../ioniz_species/Novi_thrutime/'+labels[j]+'_Novi_3456.np')
    R = np.loadtxt('../ioniz_species/Rbins_thrutime/'+labels[j]+'_Rbins_3456.np')
    plt.plot(R,Novi,label=labels[j],color=colors[j])

#plt.show()
#quit()

# Include GM4 no BHs
label_GM4_noBH = 'GM4_noBH'
color_GM4_noBH = 'Green'

f = pynbody.load('/nobackup/nnsanche/pioneer50h243GM4.1536gst1bwK1_laststep/pioneer50h243GM4.1536gst1bwK1.003456')
pynbody.analysis.halo.center(f.star)
f.physical_units()
h = f.halos()
h1 = h[1]
pynbody.analysis.angmom.faceon(h1)

r_max = 10  # kpc
twenty_kpc_incm = 6.171*(10**22)

Rg_d = ((h1.g['x'].in_units('kpc'))**2. + (h1.g['y'].in_units('kpc'))**2. + (h1.g['z'].in_units('kpc'))**2.)**(0.5)
disk_gas_xyzmax =  (Rg_d < r_max)
disk_gas_mask = disk_gas_xyzmax #& disk_gas_zmax
disk_gas = h1.g[disk_gas_mask] #& disk_gas_zmax]
CGM_gas  = h1.g[~disk_gas_mask]
CGM_temp = np.array(CGM_gas['temp'])

CGM_gas['ovi'] = pynbody.analysis.ionfrac.calculate(CGM_gas,ion='ovi',mode='new') 
m_p = 1.6726 * 10**-24 #g
print('Total mass in CGM:', np.sum(CGM_gas['mass']))
    
CGMprofile = profile.Profile(CGM_gas,min='0.1 kpc',max='250 kpc')

np.savetxt('GM4noBHs_Novi_3456.np',np.log10((CGMprofile['mass'].in_units('g')*CGMprofile['OxMassFrac']*CGMprofile['ovi']/(16*m_p))/CGMprofile._binsize.in_units('cm**2')))
np.savetxt('GM4noBHs_T_3456.np',np.log10(CGMprofile['temp'].in_units('K')))
np.savetxt('GM4noBHs_rho_3456.np',np.log10(CGMprofile['rho'].in_units('g cm^-3')))
np.savetxt('GM4noBHs_Omass_3456.np',np.log10(CGMprofile['mass'].in_units('g')*CGMprofile['OxMassFrac']/(16*m_p)))
np.savetxt('GM4noBHs_Rbins_3456.np',CGMprofile['rbins'].in_units('kpc'))

#        sph.image(CGM_gas,qty="temp",width=500,cmap="YlOrRd")
#        plt.savefig(labels[k]+'_Tmap_'+ts[i]+'.pdf')
#        plt.show()
#        plt.close()

plt.plot(CGMprofile['rbins'].in_units('kpc'),np.log10((CGMprofile['mass'].in_units('g')*CGMprofile['OxMassFrac']*CGMprofile['ovi']/(16*m_p))/CGMprofile._binsize.in_units('cm**2')),label=label_GM4_noBH,linestyle='--')

plt.title('z = 0.17')
plt.ylabel(r'log(N$_{OVI}$) [cm$^{-2}$]')
plt.xlabel('R [kpc]')
plt.ylim(12.5,16.5)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_plusGM4noBH_Novi_R.pdf')
plt.show()
plt.close()

quit()

plt.title('z = 0.0')
plt.ylabel(r'log(N$_{ovi}$) [cm$^{-2}$]')
plt.xlabel('R [kpc]')
plt.ylim(12.5,16.5)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_Novi_R.pdf')
plt.show()
plt.close()

for j in range(len(labels)-1):
    M_ox = np.loadtxt('Omass_thrutime/'+labels[j]+'_Omass_4096.np')
    R = np.loadtxt('Rbins_thrutime/'+labels[j]+'_Rbins_4096.np')
    plt.plot(R,M_ox,label=labels[j],color=colors[j])

plt.title('z = 0.0')
plt.ylabel(r'log(M$_[O]$')
plt.xlabel('R [kpc]')
#plt.ylim(12.5,16.5)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_Omass_R.pdf')
plt.show()
plt.close()

# Calculate cooling time
m_H = 1.67 * 10**-24 #g
C_1 = 3.88 * 10**11 # s K^-1/2 cm^-3
C_2 = 5 * 10**7 #K
f_m = 1.0 # metallicity dependent constant is 1 for solar metallicity Balogh et al 2013
mu = 0.6  # solar metallicity

for j in range(len(labels)-1):
    rho = np.loadtxt('rho_thrutime/'+labels[j]+'_rho_4096.np')
    R = np.loadtxt('Rbins_thrutime/'+labels[j]+'_Rbins_4096.np')
    T = np.loadtxt('T_thrutime/'+labels[k]+'_T_4096.np')
    t_cool = (C_1*mu*m_H*(10**T)**(0.5))/((10**rho)*(1+(C_2*f_m/(10**T)))) / (3.154 * 10**7)# years
    plt.plot(R,np.log10(t_cool),label=labels[j],color=colors[j])

plt.title('z = 0.0')
plt.ylabel(r't$_{cool}$ [years]')
plt.xlabel('R [kpc]')
#plt.ylim(-0.05,1.2)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_tcool_R.pdf')
plt.show()
