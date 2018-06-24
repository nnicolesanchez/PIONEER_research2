# This script confirms the ID for the most SMBH and 
# makes the following plots using dat made from
# ./BH_makefile.py: (because plotting was weird in script?)
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

for i in range(len(sims)):
    print(sims[i])
    # Read in numpy array 
    BH_times = np.loadtxt(labels[i]+'_BHtime.txt')
    BH_mass  = np.loadtxt(labels[i]+'_BHmass.txt')
    BH_mdot  = np.loadtxt(labels[i]+'_BHmdot.txt')
    
    plt.plot(BH_times,BH_mass,color=colors[i],label=labels[i])

plt.ylabel(r'M$_{BH}$/M$_{\odot}$')
plt.xlabel('t [Gyr]')
plt.show()



                        
