import pynbody.plot.sph as sph
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import os.path
import seaborn as sns
import pandas as pd
import numpy as np
import pynbody
from pynbody.analysis import profile


# Just using k = 1 and k = 2, for GM1 & GM4 for now
k = 0
## MOVED FILES FROM FABIO TO ALYSON BROOKS: /nobackupp8/ambrook2/fgoverna_pleiades_p8_files
sim = ['/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM4.1536gst1bwK1_laststep/pioneer50h243GM4.1536gst1bwK1.00']
labels = ['P0','GM1','GM4','GM5','GM6','GM7','GM4_noBH']
colors = sns.cubehelix_palette(8)
#print('LOADING SIM:',labels[k])

ts = np.loadtxt('../'+labels[k]+'/timesteps.txt',dtype=str)
time = '3456'

#for i in range(len(ts)-7,len(ts)):
for i in range(len(sim)):
    print('LOADING SIM:',labels[i],'AT TIMESTEP:',time)
#    if k == 5:
#        if ts[i] == '4096':
#            continue
    ######################
    # READ IN SIMULATION #
    ######################
    f = pynbody.load(sim[i]+time)
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

    disk_profile = profile.Profile(disk_gas,min='0.1 kpc',max='10 kpc')

    np.savetxt(labels[i]+'_metals_'+time+'.np',disk_profile['metals'])
    np.savetxt(labels[i]+'_Rbins_'+time+'.np',disk_profile['rbins'].in_units('kpc')) 

    plt.plot(disk_profile['rbins'].in_units('kpc'),disk_profile['metals'],label=labels[i],color=colors[i])

plt.title('z = 0.17')
plt.ylabel(r'$Z$')
plt.xlabel('R [kpc]')
#plt.ylim(5.2,6.2)
plt.xlim(-1,11)
plt.legend()
plt.savefig('ALLGMs_DISKmetals_R_plusGM4noBH.pdf')
plt.show()
plt.close()


quit()
#    CGM_gas  = h1.g[~disk_gas_mask]
#    CGM_temp = np.array(CGM_gas['temp'])

    #########################
    # CALCULATE OVI DENSITY #
    #########################
    # Ionization fraction of OVI compare to total Oxygen
 #   CGM_gas['ovi'] = pynbody.analysis.ionfrac.calculate(CGM_gas,ion='ovi',mode='new') 
#    m_p = 1.6726 * 10**-24 #g
#    print('Total mass in CGM:', np.sum(CGM_gas['mass']))
#    if np.sum(CGM_gas['mass']) < 1:
#        continue
#    else :

#        CGMprofile = profile.Profile(CGM_gas,min='0.1 kpc',max='250 kpc')

#        np.savetxt('Novi_thrutime/'+labels[k]+'_Novi_'+ts[i]+'.np',np.log10((CGMprofile['mass'].in_units('g')*CGMprofile['OxMassFrac']*CGMprofile['ovi']/(16*m_p))/CGMprofile._binsize.in_units('cm**2')))
#        np.savetxt('T_thrutime/'+labels[k]+'_T_'+ts[i]+'.np',np.log10(CGMprofile['temp'].in_units('K')))
#        np.savetxt('rho_thrutime/'+labels[k]+'_rho_'+ts[i]+'.np',np.log10(CGMprofile['rho'].in_units('g cm^-3')))
#        np.savetxt('Omass_thrutime/'+labels[k]+'_Omass_'+ts[i]+'.np',np.log10(CGMprofile['mass'].in_units('g')*CGMprofile['OxMassFrac']/(16*m_p)))
#        np.savetxt('Rbins_thrutime/'+labels[k]+'_Rbins_'+ts[i]+'.np',CGMprofile['rbins'].in_units('kpc'))
#        np.savetxt('totgasmass_thrutime/'+labels[k]+'_totmassgas_'+ts[i]+'.np',CGMprofile['mass'].in_units('g'))
#        np.savetxt('metals_thrutime/'+labels[k]+'_metals_'+ts[i]+'.np',CGMprofile['metals'])
        #print(CGMprofile['coolontime'].in_units('s').units)
#        np.savetxt('coolontime_thrutime/'+labels[k]+'_coolontime_'+ts[i]+'.np',CGMprofile['coolontime'].in_units('s'))
   
#        sph.image(CGM_gas,qty="temp",width=500,cmap="YlOrRd")
#        plt.savefig(labels[k]+'_Tmap_'+ts[i]+'.pdf')
#        plt.show()
#        plt.close()

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
