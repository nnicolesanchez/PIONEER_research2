# This script traces the metals in the PIONEER GM simulations.
# It creates the outputs for metal mass as a function of time for:
#      - metal mass within 1 kpc
#      - metal mass between 1 kpc & 20 kpc
#      - metal mass between 20 kpc amd the virial radius
#
# It check if these output have been created and 
# creates the plots of metal flux vs time.


# N. Nicole Sanchez -- Aug 20 2018
# Univ. of Wash.    -- Nbody Shop
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#import matplotlib as mpl
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

file_path = 'metalflux_data/'+name+'totmass_in1kpc.txt'
if (os.path.exists(file_path) == False) or ('overwrite' in sys.argv):
    print('Beginning to write data to metalflux_data/')
    steps = np.loadtxt('../'+str(sys.argv[1])+'/timesteps.txt',dtype='str')

    times = []
    totmass_in1kpc    = []
    Zmass_in1kpc      = []
    totmass_bw1_20kpc = []
    Zmass_bw1_20kpc   = []
    totmass_out20kpc  = []
    Zmass_out20kpc    = []

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
        one_kpc    = 1  # kpc
        twenty_kpc = 20
        
        h1_gas_in1kpc    = h1.g[h1.g['r'] < one_kpc]
        h1_gas_bw1_20kpc = h1.g[(h1.g['r'] > one_kpc) & (h1.g['r'] < twenty_kpc)]
        h1_gas_out20kpc  = h1.g[h1.g['r'] > twenty_kpc]

        print('Total gas mass < 1 kpc:', np.sum(h1_gas_in1kpc['mass'].in_units("Msol")))
        print('Total metal mass < 1 kpc:',np.sum(h1_gas_in1kpc['mass'].in_units("Msol")*h1_gas_in1kpc['metals']))
        totmass_in1kpc.append(np.sum(h1_gas_in1kpc['mass'].in_units("Msol")))
        Zmass_in1kpc.append(np.sum(h1_gas_in1kpc['mass'].in_units("Msol")*h1_gas_in1kpc['metals']))
        
        print('Total gas mass > 1 kpc & < 20 kpc:', np.sum(h1_gas_bw1_20kpc['mass'].in_units("Msol")))
        print('Total metal mass > 1 kpc & < 20 kpc:',np.sum(h1_gas_bw1_20kpc['mass'].in_units("Msol")*h1_gas_bw1_20kpc['metals'])) 
        totmass_bw1_20kpc.append(np.sum(h1_gas_bw1_20kpc['mass'].in_units("Msol")))
        Zmass_bw1_20kpc.append(np.sum(h1_gas_bw1_20kpc['mass'].in_units("Msol")*h1_gas_bw1_20kpc['metals']))
        
        print('Total gas mass > 20 kpc:', np.sum(h1_gas_out20kpc['mass'].in_units("Msol")))
        print('Total metal mass > 20 kpc:',np.sum(h1_gas_out20kpc['mass'].in_units("Msol")*h1_gas_out20kpc['metals'])) 
        totmass_out20kpc.append(np.sum(h1_gas_out20kpc['mass'].in_units("Msol")))
        Zmass_out20kpc.append(np.sum(h1_gas_out20kpc['mass'].in_units("Msol")*h1_gas_out20kpc['metals']))
        
        times.append(sim.properties['time'].in_units("Gyr"))
        
        np.savetxt('metalflux_data/'+name+'totmass_in1kpc.txt',totmass_in1kpc)
        np.savetxt('metalflux_data/'+name+'Zmass_in1kpc.txt',Zmass_in1kpc)
        np.savetxt('metalflux_data/'+name+'totmass_bw1_20kpc.txt',totmass_bw1_20kpc)
        np.savetxt('metalflux_data/'+name+'Zmass_bw1_20kpc.txt',Zmass_bw1_20kpc)
        np.savetxt('metalflux_data/'+name+'totmass_out20kpc.txt',totmass_out20kpc)
        np.savetxt('metalflux_data/'+name+'Zmass_out20kpc.txt',Zmass_out20kpc)
        np.savetxt('metalflux_data/'+name+'times.txt',times)
        print(' ')

else:
    print('Reading in files from metalflux_data/')
    totmass_in1kpc    = np.loadtxt('metalflux_data/'+name+'totmass_in1kpc.txt',unpack=True)
    Zmass_in1kpc      = np.loadtxt('metalflux_data/'+name+'Zmass_in1kpc.txt',unpack=True)
    totmass_bw1_20kpc = np.loadtxt('metalflux_data/'+name+'totmass_bw1_20kpc.txt',unpack=True)
    Zmass_bw1_20kpc   = np.loadtxt('metalflux_data/'+name+'Zmass_bw1_20kpc.txt',unpack=True)
    totmass_out20kpc  = np.loadtxt('metalflux_data/'+name+'totmass_out20kpc.txt',unpack=True)
    Zmass_out20kpc    = np.loadtxt('metalflux_data/'+name+'Zmass_out20kpc.txt',unpack=True)
    times = np.loadtxt('metalflux_data/'+name+'times.txt',unpack=True)

# Individual plots
#fig = plt.figure(figsize=(10, 12))
#ax1 = fig.add_subplot(311)
#ax2 = fig.add_subplot(312)
#ax3 = fig.add_subplot(313)
#ax1.plot(times,np.log10(Zmass_in1kpc))
#ax2.plot(times,np.log10(Zmass_bw1_20kpc))
#ax3.plot(times,np.log10(Zmass_out20kpc))
#ax1.set_ylabel(r'Z Mass < 1 kpc [M$_{\odot}$]')
#ax2.set_ylabel(r'Z Mass > 1 kpc & < 20 kpc [M$_{\odot}$]')
#ax3.set_ylabel(r'Z Mass > 20 kpc [M$_{\odot}$]')
#ax3.set_xlabel('Time [Gyr]')
#plt.show()

plt.plot(times,np.log10(Zmass_in1kpc),label=r'r < 1 kpc',linestyle=':',color='Pink')
plt.plot(times,np.log10(Zmass_bw1_20kpc),label=r'1 kpc < r <  20 kpc',linestyle='--',color='Salmon')
plt.plot(times,np.log10(Zmass_out20kpc),label=r'20 kpc < r < R$_{vir}$',linestyle='-',color='SteelBlue')
plt.ylim(3,9)
plt.xlabel('Age [Gyr]')
plt.ylabel(r'log M$_Z$ [M$_{\odot}$]')
plt.legend()
plt.savefig(name+'_Zmass_thrutime.pdf')
plt.show()
