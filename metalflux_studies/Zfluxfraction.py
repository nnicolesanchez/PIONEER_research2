# This script traces the metal FRACTIONS in the PIONEER GM simulations.
# It creates the outputs for metal mass FRACTION as a function of time for:
#      - metal mass within 10 kpc
#      - metal mass between 10 kpc amd the virial radius
#
# It check if these output have been created and 
# creates the plots of:
#      (total metal within 10 kpc) - (total metal mass outside 10 kpc)/
#                    total metal mass in timestep
# Plots this as a function of time.

# N. Nicole Sanchez -- Aug 23 2018
# Univ. of Wash.    -- Nbody Shop
import matplotlib.pyplot as plt
import numpy as np
import pynbody
import os.path
import sys


ALL_sims = ['/nobackupp2/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM1.1536gst1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00','/nobackupp2/nnsanche/NO_BHs/pioneer50h243.1536gst1bwK1/pioneer50h243.1536gst1bwK1.00','/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM1.1536gst1bwK1/pioneer50h243GM1.1536gst1bwK1.00','/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM7.1536gst1bwK1/pioneer50h243GM7.1536gst1bwK1.00','/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM4.1536gst1bwK1/pioneer50h243GM4.1536gst1bwK1.00']

if len(sys.argv) == 1:
    print('No galaxy selected. Current options: P0, GM1, GM7 (aka GM2), GM4 (aka GM3)')
    print('Includes galaxies without BH physics: Call "P0noBH"')
    print('Syntax: "Zfluxfraction.py GM1noBH"')
    quit()
else:
    if (str(sys.argv[1]) == 'P0'):
        k = 0
    elif (str(sys.argv[1]) == 'GM1'):
        k = 1
    elif (str(sys.argv[1]) == 'GM7'):
        k = 2
    elif (str(sys.argv[1]) == 'GM4'):
        k = 3
    elif (str(sys.argv[1]) == 'P0noBH'):
        k = 4
    elif (str(sys.argv[1]) == 'GM1noBH'):
        k = 5
    elif (str(sys.argv[1]) == 'GM7noBH'):
        k = 6
    elif (str(sys.argv[1]) == 'GM4noBH'):
        k = 7

    else :
        print('Not a valid option. Current options: P0, GM1, GM7 (aka GM2), GM4 (aka GM3)')
        print('Includes galaxies without BH physics: Call "P0noBH"')
        print('Syntax: "metalflux_thrutime.py GM1noBH"')
        quit()    

    if str(sys.argv[1]) == 'GM7':
        name = 'GM2'
    elif str(sys.argv[1]) == 'GM4':
        name = 'GM3'
    elif str(sys.argv[1]) == 'GM7noBH':
        name = 'GM2noBH'
    elif str(sys.argv[1]) == 'GM4noBH':
        name = 'GM3noBH'
    else:
        name = str(sys.argv[1])

file_path = 'metalflux_data/'+name+'_Zflux_times.py'
if (os.path.exists(file_path) == False) or ('overwrite' in sys.argv):
    print('Beginning to write data to metalflux_data/')
    steps = np.loadtxt('../'+str(sys.argv[1])+'/timesteps.txt',dtype='str')

    times = []
    Zmass_inner = []
    Zmass_outer = []
    totZmass    = []

    for ts in range(0,len(steps)):
        sim   = pynbody.load(ALL_sims[k]+steps[ts])
        print(name+' simulation at z = ','%.2f' % sim.properties['z'],' and time = ',sim.properties['time'].in_units("Gyr"))

        #sim.loadable_keys()
        sim.physical_units()
        h = sim.halos()
        h1 = h[1]
        pynbody.analysis.angmom.faceon(h1)
        
        # Constants
        m_H = 1.6733 * 10**-24 #g
        Z_sun = 0.0142 # (Asplund 2009; https://arxiv.org/pdf/0909.0948.pdf) 
        
        # Want to isolate metals in:
        # Isolate gas within 1 kpc, between 1-20 kpc, and outside 20 kpc to Rvir
        bound = 10 #20 #kpc
        
        h1_gas_inner = h1.g[h1.g['r'] < bound]
        h1_gas_outer = h1.g[h1.g['r'] >= bound]
        
        print('Total gas mass < 10 kpc:', np.sum(h1_gas_inner['mass'].in_units("Msol")))
        print('Total metal mass < 10 kpc:',np.sum(h1_gas_inner['mass'].in_units("Msol")*h1_gas_inner['metals'])) 
        totmass_inner.append(np.sum(h1_gas_inner['mass'].in_units("Msol")))
        Zmass_inner.append(np.sum(h1_gas_inner['mass'].in_units("Msol")*h1_gas_inner['metals']))
        
        print('Total gas mass > 20 kpc:', np.sum(h1_gas_outer['mass'].in_units("Msol")))
        print('Total metal mass > 20 kpc:',np.sum(h1_gas_outer['mass'].in_units("Msol")*h1_gas_outer['metals'])) 
        totmass_outer.append(np.sum(h1_gas_outer['mass'].in_units("Msol")))
        Zmass_outer.append(np.sum(h1_gas_outer['mass'].in_units("Msol")*h1_gas_outer['metals']))
        
        times.append(sim.properties['time'].in_units("Gyr"))
        
        np.savetxt('metalflux_data/'+name+'totmass_inner.txt',totmass_inner)
        np.savetxt('metalflux_data/'+name+'Zmass_inner.txt',Zmass_inner)
        np.savetxt('metalflux_data/'+name+'totmass_outer.txt',totmass_outer)
        np.savetxt('metalflux_data/'+name+'Zmass_outer.txt',Zmass_outer)
        np.savetxt('metalflux_data/'+name+'times.txt',times)
        print(' ')

else:
    print('Reading in files from metalflux_data/')
    totmass_inner = np.loadtxt('metalflux_data/'+name+'totmass_inner.txt',unpack=True)
    Zmass_inner   = np.loadtxt('metalflux_data/'+name+'Zmass_inner.txt',unpack=True)
    totmass_outer  = np.loadtxt('metalflux_data/'+name+'totmass_outer.txt',unpack=True)
    Zmass_outer    = np.loadtxt('metalflux_data/'+name+'Zmass_outer.txt',unpack=True)
    times = np.loadtxt('metalflux_data/'+name+'times.txt',unpack=True)


plt.plot(BH_data.calculate_for_progenitors('t()')[0][:m],np.log10(BH_data.calculate_for_progenitors('BH_mass')[0]),color='DodgerBlue',label='BH mass',linewidth=2)
plt.plot(times,np.log10(Zmass_inner),label=r'Mass of gas in metals within 20 kpc',linestyle='--',color='Salmon')
plt.plot(times,np.log10(Zmass_outer),label=r'Mass of gas in metals with 20 kpc < r < R$_{vir}$',linestyle=':',color='SteelBlue')
plt.ylim(6,9)
plt.xlabel('Age/Gyr',fontsize=15)
plt.ylabel(r'log M/M$_{\odot}$',fontsize=15)
plt.text(0.2,8.8,name,size=15)
plt.legend()
plt.savefig(name+'_Zmass_BHmass_thrutime.pdf')
plt.show()




