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

totmetals_array = []
CGMmetals_array = []

ts       = np.loadtxt('../'+labels[k]+'/timesteps.txt',dtype=str)
time_Gyr = np.loadtxt('../'+labels[k]+'/times_gyr.txt',dtype=str)  
redshift = np.loadtxt('../'+labels[k]+'/redshifts.txt',dtype=str) 
for i in range(5,len(ts)):
    print('TIMESTEP:',ts[i],'TIME:',time_Gyr[i],'z =',redshift[i])
    if k == 5:
        if ts[i] == '4096':
            continue

    Z     = np.loadtxt('TOTGXY_metals_thrutime/'+labels[k]+'_metals_'+ts[i]+'.np')
    Rbins = np.loadtxt('TOTGXY_Rbins_thrutime/'+labels[k]+'_Rbins_'+ts[i]+'.np')
    plt.plot(Rbins,Z)
    plt.title('z = '+redshift[i])
    plt.ylabel(r'$Z/Z_{\odot}$')
    plt.xlabel('R [kpc]')
    plt.ylim(0,0.05)
    plt.xlim(0,250)
    plt.savefig(labels[k]+'_plots/TOTGXY_metals_Rbins_'+ts[i]+'.pdf')
    plt.clf()
#    plt.show()

    totmetals_array.append(np.mean(Z))

np.savetxt(labels[k]+'_meanmetals_time.txt',totmetals_array)

plt.plot(time_Gyr[5:],totmetals_array)
plt.xlabel('Time [Gyr]')
plt.ylabel(r'$Z/Z_{\odot}$')
plt.savefig(labels[k]+'_meanmetals_time.pdf')
plt.show()

#    CGMmetals_array.append(np.mean(CGMprofile['metals']))

#    np.savetxt('CGM_Novi_thrutime/'+labels[k]+'_Novi_'+ts[i]+'.np',np.log10((CGMprofile['mass'].in_units('g')*CGMprofile['OxMassFrac']*CGMprofile['ovi']/(16*m_p))/CGMprofile._binsize.in_units('cm**2')))
#    np.savetxt('CGM_T_thrutime/'+labels[k]+'_T_'+ts[i]+'.np',np.log10(CGMprofile['temp'].in_units('K')))
#    np.savetxt('CGM_rho_thrutime/'+labels[k]+'_rho_'+ts[i]+'.np',np.log10(CGMprofile['rho'].in_units('g cm^-3')))
#    np.savetxt('CGM_Omass_thrutime/'+labels[k]+'_Omass_'+ts[i]+'.np',np.log10(CGMprofile['mass'].in_units('g')*CGMprofile['OxMassFrac']/(16*m_p)))
#    np.savetxt('CGM_Rbins_thrutime/'+labels[k]+'_Rbins_'+ts[i]+'.np',CGMprofile['rbins'].in_units('kpc'))
#    np.savetxt('CGM_totgasmass_thrutime/'+labels[k]+'_totmassgas_'+ts[i]+'.np',CGMprofile['mass'].in_units('g'))
#    np.savetxt('CGM_metals_thrutime/'+labels[k]+'_metals_'+ts[i]+'.np',CGMprofile['metals'])
   
