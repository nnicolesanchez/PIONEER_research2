# This script plots the star formation history for the GM
# suite for h243 (P0, GM1, GM2(GM7), GM3(GM4)
# With and without black hole physics

# N. Nicole Sanchez  -- Created: Spring 2017
# Univ of W, Seattle -- Edited: June 24, 2018
import matplotlib.pyplot as plt
import numpy as np
import pynbody
import os

# Information for loading all available sims P0 - GM7
sims      = ['/nobackupp2/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.004096','/nobackupp2/nnsanche/pioneer50h243GM1.1536gst1bwK1BH_Simpyorbitnotworking/pioneer50h243GM1.1536gst1bwK1BH.004096','/nobackupp2/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.003968','/nobackupp2/nnsanche/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.004096','/nobackupp2/nnsanche/NO_BHs/pioneer50h243.1536gst1bwK1/pioneer50h243.1536gst1bwK1.003456','/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM1.1536gst1bwK1/pioneer50h243GM1.1536gst1bwK1.003456','/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM7.1536gst1bwK1/pioneer50h243GM7.1536gst1bwK1.003456','/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM4.1536gst1bwK1/pioneer50h243GM4.1536gst1bwK1.003456']

labels = ['P0','GM1','GM2','GM3','P0_noBH','GM1_noBH','GM2_noBH','GM3_noBH']
colors = ['DodgerBlue','SteelBlue','FireBrick','Salmon','DodgerBlue','SteelBlue','FireBrick','Salmon']
lines  = ['-','-','-','-','--','--','--','--']

# Create sfh data arrays (also save 
for i in range(0,4):#,len(sims)):
    print(labels[i])

    if os.path.exists(labels[i]+'_sfhistory_bins.txt'):
        sfhist, bins = np.genfromtxt(labels[i]+'_sfhistory_bins.txt',unpack=True)
        plt.plot(bins,sfhist,color=colors[i],linestyle=lines[i],label=labels[i],drawstyle='steps')
    else:
        sim = pynbody.load(sims[i])
        h   = sim.halos()
        h1  = h[1]
        pynbody.analysis.halo.center(h1,mode='ssc')
        pynbody.analysis.angmom.faceon(h1)
        sim.physical_units()

        # Get sfh info for whole sim; change "sim" below to "h1" for Main Halo sfh
        sfhist, bins = pynbody.plot.stars.sfh(sim,massform=True,color=colors[i],linestyle=lines[i],legend=True,trange=[0,14],label=labels[i])

        # Save for easier plot making later (tho this script doesn't take that long to run)
        # Good practice to keep in place though
        np.savetxt(labels[i]+'_sfhistory_bins.txt', np.transpose([sfhist, bins[:-1]]))

        plt.plot(bins[:-1],sfhist,color=colors[i],linestyle=lines[i],label=labels[i])

# For individual labels instead of legend (bc Fabio asked. Shrug. I guess it looks better?)
# Remember to "unset" legend when using these instead
#plt.text(10,11,'P0',color=colors[0])
#plt.text(8,7.5,'GM1',color=colors[1])
#plt.text(4,1,'GM4',color=colors[2])
#plt.text(3.5,6,'GM5',color=colors[3])
#plt.text(2,4,'GM6',color=colors[4])
#plt.text(1,3,'GM7',color=colors[5])

plt.ylabel(r'SFR [M$_{\odot}$ yr$^{-1}$]',size=15)
plt.xlabel('Time [Gyr]',size=15)
#plt.ylim(0,30) #noBH
plt.ylim(0,18) #BH
plt.xlim(0,14)
plt.legend(fontsize=15)
plt.savefig('ALLBH_sfh.pdf')
plt.show()
