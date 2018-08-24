# This script pulls out the following properties for the 
# ChaNGa galaxies: Patient 0, GM1, GM4, GM5, GM6
#        - Stellar mass of galaxy
#        - Gas mass of galaxy
#        * Split into CGM vs Disk?
#        - Merger? yes no
#            - Mass ratio if yes
#            - Redshift if yes

# N. Nicole Sanchez -- June 29 2017
# Univ. of Wash.    -- Edited: June 24, 2018
import matplotlib.pyplot as plt
import numpy as np
import pynbody
import sys

timesteps = ['1280','1408','1536','1664','1739'] # two before and two after merger!
sims = ['/nobackupp2/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM1.1536gst1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.00']
labels = ['P0','GM1','GM4','GM5','GM6','GM7']

for i in range(len(sims)-2):
    #timesteps = np.loadtxt('../'+labels[i]+'/timesteps.txt',dtype='str')
    if i == 0 or i == 3:
        continue
    else:
        for j in range(len(timesteps)):
            
            print('Making plot for',labels[i],' at ',timesteps[j])
            sim = pynbody.load(sims[i]+timesteps[j],verbose=True)
            sim.physical_units()
            h = sim.halos()
            h1 = h[1]
            pynbody.analysis.angmom.faceon(h1)
            
            pynbody.plot.sph.image(h1.g,qty='temp',width='20 kpc',vmin=0.1, vmax=10**6)
            plt.title(labels[i]+', z ='+timesteps[j])
            plt.savefig(labels[i]+'_xy_temp_20kpc_'+timesteps[j]+'_grp.pdf')
#            plt.show()
    

    

        

