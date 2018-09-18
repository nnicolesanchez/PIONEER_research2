# This script plots the star formation history for the GM
# suite for h243 (all on one plot)
import matplotlib.pyplot as plt
import numpy as np
import pynbody



sims_noBHs = ['/nobackup/nnsanche/NO_BHs/pioneer50h243.1536gst1bwK1/pioneer50h243.1536gst1bwK1.004096','/nobackup/nnsanche/NO_BHs/pioneer50h243GM1.1536gst1bwK1_3456/pioneer50h243GM1.1536gst1bwK1.003456','/nobackup/nnsanche/NO_BHs/pioneer50h243GM4.1536gst1bwK1/pioneer50h243GM4.1536gst1bwK1.004096','/nobackup/nnsanche/NO_BHs/pioneer50h243GM5.1536gst1bwK1/pioneer50h243GM5.1536gst1bwK1.003456','/nobackup/nnsanche/NO_BHs/pioneer50h243GM6.1536gst1bwK1/pioneer50h243GM6.1536gst1bwK1.003712','/nobackup/nnsanche/NO_BHs/pioneer50h243GM7.1536gst1bwK1/pioneer50h243GM7.1536gst1bwK1.003456']
# Include GM4 no BHs
labels_noBHs = ['P0_noBH','GM1_noBH','GM4_noBH','GM5_noBH','GM6_noBH','GM7_noBH']
colors_noBHs = ['DodgerBlue','SteelBlue','FireBrick','IndianRed','Salmon','Orange']
linestyle = ['--']



# Create sfh data arrays (also save 
for i in range(len(sims_noBHs)):
    print(labels_noBHs[i])
    sim = pynbody.load(sims_noBHs[i])
    h = sim.halos()
    h1 = h[1]
    pynbody.analysis.halo.center(h1,mode='ssc')
    pynbody.analysis.angmom.faceon(h1)
    sim.physical_units()

    print('Halo mass:',h1['mass'].sum())
    # Get sfh info for whole sim; change "sim" below to "h1" for Main Halo sfh
    sfhist, bins = pynbody.plot.stars.sfh(sim,filename=labels[i]+'sfh_xlim.pdf',massform=True,color=colors[i],legend=True,trange=[0,14],label=labels[i])


    # Save for easier plot making later (tho this script doesn't take that long to run)
    # Good practice to keep in place though
    np.savetxt(labels_noBHs[i]+'_sfhistory_bins.txt', np.transpose([sfhist, bins[:-1]]))

# For individual labels instead of legend (bc Fabio asked. Shrug. I guess it looks better?)
# Remember to "unset" legend when using these instead
#plt.text(10,11,'P0',color=colors[0])
#plt.text(8,7.5,'GM1',color=colors[1])
#plt.text(4,1,'GM4',color=colors[2])
#plt.text(3.5,6,'GM5',color=colors[3])
#plt.text(2,4,'GM6',color=colors[4])
#plt.text(1,3,'GM7',color=colors[5])

plt.title('z = 0.17')
plt.ylim(0,25)
#plt.xlim(5,7)
plt.legend(loc=2)
plt.savefig('Test_newnoBH.pdf')
#plt.savefig('P0-GM7_sfh_plusGM4noBH.pdf')
plt.show()

