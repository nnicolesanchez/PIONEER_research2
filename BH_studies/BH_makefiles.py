# This script confirms the ID for the most SMBH and 
# makes the following plots using Simpy:
#       - BH mass vs time
#       - BH accretion rate

# N Nicole Sanchez -- Created: June 24, 2018
# Univ. W, Seattle -- Edited: June 24, 2018
import matplotlib.pyplot as plt
import numpy as np
import Simpy


# Information for loading all available sims P0 - GM7
sims  = ['/nobackupp2/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH','/nobackupp2/nnsanche/pioneer50h243GM1.1536gst1bwK1BH_alysonver/pioneer50h243GM1.1536gst1bwK1BH','/nobackupp2/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH','/nobackupp2/nnsanche/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH']

labels = ['P0','GM1','GM2','GM3']
colors = ['DodgerBlue','SteelBlue','FireBrick','Salmon']
lines  = ['-','-','-','-']

i = 1
for i in range(len(sims)):
    print(sims[i])
    # Create the BHorbit object
    bhorbit = Simpy.BlackHoles.orbit.Orbit(sims[i], savefile=labels[i]+'_bhorbit.pkl')
    print(bhorbit.data)
    print('BHs in sim:',len(bhorbit.bhiords))
    print('BHs iords:',bhorbit.bhiords)
    
    # All the BHs in the simulation
    bh_iords = np.array(bhorbit.bhiords)

    # Sanity check to make sure index 0 is the most massive BH: It is.
    #for bh in range(len(bh_iords)):
    #    print(bhorbit.single_BH_data(bh_iords[bh], 'mass')[-1])
    
    BH_times = np.array(bhorbit.single_BH_data(bh_iords[0], 'time'))
    BH_mass  = np.array(bhorbit.single_BH_data(bh_iords[0], 'mass'))
    BH_mdot  = np.array(bhorbit.single_BH_data(bh_iords[0], 'mdot'))
    print(BH_times,BH_mass,BH_mdot)
    
    np.savetxt(labels[i]+'_BHtime.txt',BH_times)
    np.savetxt(labels[i]+'_BHmass.txt',BH_mass)
    np.savetxt(labels[i]+'_BHmdot.txt',BH_mdot)
    

                        

