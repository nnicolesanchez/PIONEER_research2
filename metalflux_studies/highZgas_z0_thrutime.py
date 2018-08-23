# This script traces the metals in the PIONEER GM simulations.
# It find the particles with the highest metallicity at z = 0
# and traces them backwards through time to their initial r.
#      - writes R (distance from center for each particle at all z)
#
# It check if these output have been created and 
# creates histogram of distances at each z


# N. Nicole Sanchez -- Aug 21 2018
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
    print('Syntax: "metalflux_thrutime.py GM1noBH"')
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

file_path = 'metalflux_data/'+name+'_highZ_r_4096.txt'
if (os.path.exists(file_path) == False) or ('overwrite' in sys.argv):
    print('Beginning to write data to metalflux_data/')
    steps = np.loadtxt('../'+str(sys.argv[1])+'/timesteps.txt',dtype='str')
    steps = steps[::-1]

    ts = 0
    sim   = pynbody.load(ALL_sims[k]+steps[ts])
    print(name+' simulation at z = ','%.2f' % sim.properties['z'],' and time = ',sim.properties['time'].in_units("Gyr"))
    
    #sim.loadable_keys()
    sim.physical_units()
    h = sim.halos()
    h1 = h[1]
    pynbody.analysis.angmom.faceon(h1)
    Rvir = float(pynbody.analysis.halo.virial_radius(sim))
    
    # Constants
    m_H = 1.6733 * 10**-24 #g
    Z_sun = 0.0142 # (Asplund 2009; https://arxiv.org/pdf/0909.0948.pdf) 
    
    # Want to isolate high metals in CGM gas:
    # Isolate gas within 1 kpc, between 1-20 kpc, and outside 20 kpc to Rvir
    ten_kpc = 10 #kpc
    
    CGM_gas              = h1.g[h1.g['r'] > ten_kpc]
    CGM_gas['ZoverZsun'] = CGM_gas['metals']/Z_sun
    highZ_CGM_gas        = CGM_gas[CGM_gas['ZoverZsun'] > 1]
    
    print('Total gas mass in CGM',np.sum(CGM_gas['mass'].in_units("Msol")))
    print('Total gas mass with Z > 1 Zsun:', np.sum(highZ_CGM_gas['mass'].in_units("Msol")))
    print('Number of gas particles with Z > 1 Zsun',len(highZ_CGM_gas['mass']))
          
    pynbody.plot.sph.image(highZ_CGM_gas,qty='rho',width=str(Rvir)+" kpc")        
    plt.savefig('plots_highZ_thrutime/'+name+'_highZ_rho_'+steps[ts]+'.pdf')
    plt.title(name+', '+ "%.1f" % sim.properties['time'].in_units("Gyr")+' Gyr')
#    plt.show()

    np.savetxt('metalflux_data/'+name+'_highZ_iord_'+steps[ts]+'.txt',highZ_CGM_gas['iord'])
    np.savetxt('metalflux_data/'+name+'_highZ_mass_'+steps[ts]+'.txt',highZ_CGM_gas['mass'])
    np.savetxt('metalflux_data/'+name+'_highZ_r_'+steps[ts]+'.txt',highZ_CGM_gas['r'])

    plt.hist(highZ_CGM_gas['r'],bins=50,log=True)
    plt.title(name+', '+ "%.1f" % sim.properties['time'].in_units("Gyr"))
    plt.xlabel('Distance from Galactic Center [kpc]')
    plt.ylabel(r'dN/dr')
    plt.ylim(0,22000)
    plt.savefig('plots_highZ_thrutime/hist_'+name+'_Zmass_'+steps[ts]+'.pdf')
#    plt.show()

file_path = 'metalflux_data/'+name+'_highZ_r_1739.txt'
if (os.path.exists(file_path) == False) or ('overwrite' in sys.argv):
    steps = np.loadtxt('../'+str(sys.argv[1])+'/timesteps.txt',dtype='str')
    steps = steps[::-1]

    # Now that we have found the highZ gas at 
    for ts in range(1,len(steps)):
        sim   = pynbody.load(ALL_sims[k]+steps[ts])
        print(name+' simulation at z = ','%.2f' % sim.properties['z'],' and time = ',sim.properties['time'].in_units("Gyr"))

        sim.physical_units()
        h = sim.halos()
        h1 = h[1]
        pynbody.analysis.angmom.faceon(h1)

        
        highZ_CGM_gas        = CGM_gas[CGM_gas['iord']]
        
        #print('Total gas mass in CGM',np.sum(CGM_gas['mass'].in_units("Msol")))
        print('Total gas mass with Z > 1 Zsun:', np.sum(highZ_CGM_gas['mass'].in_units("Msol")))
        print('Number of gas particles with Z > 1 Zsun',len(highZ_CGM_gas['mass']))
        
        pynbody.plot.sph.image(highZ_CGM_gas,qty='rho',width=str(Rvir)+" kpc")
        plt.title(name+', '+ "%.1f" % sim.properties['time'].in_units("Gyr")+' Gyr')
        plt.savefig('plots_highZ_thrutime/'+name+'_highZ_rho_'+steps[ts]+'.pdf')
#        plt.show()
        
        np.savetxt('metalflux_data/'+name+'_highZ_iord_'+steps[ts]+'.txt',highZ_CGM_gas['iord'])
        np.savetxt('metalflux_data/'+name+'_highZ_mass_'+steps[ts]+'.txt',highZ_CGM_gas['mass'])
        np.savetxt('metalflux_data/'+name+'_highZ_mass_'+steps[ts]+'.txt',highZ_CGM_gas['r'])

        plt.hist(highZ_CGM_gas['r'],bins=50,log=True)
        plt.title(name+', '+ "%.1f" % sim.properties['time'].in_units("Gyr")+'Gyr')
        plt.xlabel('Distance from Galactic Center [kpc]')
        plt.ylabel(r'dN/dr')
        plt.ylim(0,22000)
        plt.savefig('plots_highZ_thrutime/hist_'+name+'_Zmass_'+steps[ts]+'.pdf')
#        plt.show()

else:
    print('Reading in files from metalflux_data/')
    steps = np.loadtxt('../'+str(sys.argv[1])+'/timesteps.txt',dtype='str')
    steps = steps[::-1]
    times = np.loadtxt('metalflux_data/'+name+'times.txt')
    
    highZ_CGM_gas = []
    for ts in range(0,len(steps)):
        highZ_CGM_gas['iord'] = np.loadtxt('metalflux_data/'+name+'_highZ_iord_'+steps[ts]+'.txt')
        highZ_CGM_gas['mass'] = np.savetxt('metalflux_data/'+name+'_highZ_mass_'+steps[ts]+'.txt')
        highZ_CGM_gas['r']    = np.savetxt('metalflux_data/'+name+'_highZ_mass_'+steps[ts]+'.txt')
        

        plt.hist(highZ_CGM_gas['r'],bins=50,log=True)
        plt.title(name+', '+ "%.1f" % times[ts]+' Gyr')
        #plt.ylim(3,9)
        plt.xlabel('Distance from Galactic Center [kpc]')
        plt.ylabel(r'dN/dr')
        plt.ylim(0,22000)
        plt.savefig('plots_highZ_thrutime/hist_'+name+'_Zmass_'+steps[ts]+'.pdf')
#        plt.show()
