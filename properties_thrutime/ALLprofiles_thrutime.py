from pynbody.analysis import profile
import pynbody.plot.sph as sph
import matplotlib.pyplot as plt
import sys
import numpy as np
import pynbody


if len(sys.argv) == 1:
    print('No galaxy selected. Current options include: P0, GM1, GM4, GM5, GM6, GM7')
    print('Syntax: "ALLprofiles_thrutime.py P0"')
    quit()
elif (str(sys.argv[1]) == 'P0'):
    k = 0
elif (str(sys.argv[1]) == 'GM1'):
    k = 1
elif (str(sys.argv[1]) == 'GM4'):
    k = 2
elif (str(sys.argv[1]) == 'GM5'):
    k = 3
elif (str(sys.argv[1]) == 'GM6'):
    k = 4
elif (str(sys.argv[1]) == 'GM7'):
    k = 5

sim = ['/nobackupp2/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.00']
labels = ['P0','GM1','GM4','GM5','GM6','GM7']
print('LOADING SIM:',labels[k])

totmetals_array   = []
CGMmetals_array   = []
time_array     = []
redshift_array = []

ts = np.loadtxt('../'+labels[k]+'/timesteps.txt',dtype=str)
for i in range(0,len(ts)):
    print('LOADING TIMESTEP:',ts[i])
    if k == 5:
        if ts[i] == '4096':
            continue
    ######################
    # READ IN SIMULATION #
    ######################
    f = pynbody.load(sim[k]+ts[i])
    print('Time:',f.properties['time'].in_units('Gyr'),'Redshift',f.properties['z'])
    pynbody.analysis.halo.center(f.star)
    f.physical_units()
    h = f.halos()
    h1 = h[1]
    pynbody.analysis.angmom.faceon(h1)
    
    if np.sum(h1.g['mass']) < 1:
        continue
    else :

        TOTGXYprofile = profile.Profile(h1.g,min='0.1 kpc',max='250 kpc')

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
    if np.sum(CGM_gas['mass']) < 1:
        continue
    else :

        CGMprofile = profile.Profile(CGM_gas,min='0.1 kpc',max='250 kpc')

        np.savetxt('TOTGXY_metals_thrutime/'+labels[k]+'_metals_'+ts[i]+'.np',TOTGXYprofile['metals'])
        np.savetxt('TOTGXY_Rbins_thrutime/'+labels[k]+'_Rbins_'+ts[i]+'.np',TOTGXYprofile['rbins'])
        np.savetxt('CGM_Novi_thrutime/'+labels[k]+'_Novi_'+ts[i]+'.np',np.log10((CGMprofile['mass'].in_units('g')*CGMprofile['OxMassFrac']*CGMprofile['ovi']/(16*m_p))/CGMprofile._binsize.in_units('cm**2')))
        np.savetxt('CGM_T_thrutime/'+labels[k]+'_T_'+ts[i]+'.np',np.log10(CGMprofile['temp'].in_units('K')))
        np.savetxt('CGM_rho_thrutime/'+labels[k]+'_rho_'+ts[i]+'.np',np.log10(CGMprofile['rho'].in_units('g cm^-3')))
        np.savetxt('CGM_Omass_thrutime/'+labels[k]+'_Omass_'+ts[i]+'.np',np.log10(CGMprofile['mass'].in_units('g')*CGMprofile['OxMassFrac']/(16*m_p)))
        np.savetxt('CGM_Rbins_thrutime/'+labels[k]+'_Rbins_'+ts[i]+'.np',CGMprofile['rbins'].in_units('kpc'))
        np.savetxt('CGM_totgasmass_thrutime/'+labels[k]+'_totmassgas_'+ts[i]+'.np',CGMprofile['mass'].in_units('g'))
        np.savetxt('CGM_metals_thrutime/'+labels[k]+'_metals_'+ts[i]+'.np',CGMprofile['metals'])
        #np.savetxt('coolontime_thrutime/'+labels[k]+'_coolontime_'+ts[i]+'.np',CGMprofile['coolontime'].in_units('s'))
   
        totmetals_array.append(np.mean(TOTGXYprofile['metals']))
        CGMmetals_array.append(np.mean(CGMprofile['metals']))
        time_array.append(f.properties['time'].in_units('Gyr'))
        redshift_array.append(f.properties['z'])
