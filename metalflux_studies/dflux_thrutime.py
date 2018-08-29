# This script reads in the data created by metalflux_thrutime.py
# Creates a plot of the metal flux through 10 kpc
#      - total mass of inner metals - total mass of outer metals
#                  / total metals in timestep

# N. Nicole Sanchez -- Aug 20 2018
# Univ. of Wash.    -- Nbody Shop 
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import numpy as np
import pynbody
import os.path
import tangos
import sys


ALL_sims = ['/nobackupp2/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM1.1536gst1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00','/nobackupp2/nnsanche/NO_BHs/pioneer50h243.1536gst1bwK1/pioneer50h243.1536gst1bwK1.00','/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM1.1536gst1bwK1/pioneer50h243GM1.1536gst1bwK1.00','/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM7.1536gst1bwK1/pioneer50h243GM7.1536gst1bwK1.00','/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM4.1536gst1bwK1/pioneer50h243GM4.1536gst1bwK1.00']
name   = ['P0','GM1','GM2','GM3']
labels = ['P0','GM1','GM2','GM3']
colors = ['SteelBlue','DodgerBlue','FireBrick','Salmon']

print('Reading in files from metalflux_data/')
for i in range(len(name)):
    totmass_inner = np.loadtxt('metalflux_data/'+name[i]+'totmass_inner.txt',unpack=True)
    Zmass_inner   = np.loadtxt('metalflux_data/'+name[i]+'Zmass_inner.txt',unpack=True)
    totmass_outer  = np.loadtxt('metalflux_data/'+name[i]+'totmass_outer.txt',unpack=True)
    Zmass_outer    = np.loadtxt('metalflux_data/'+name[i]+'Zmass_outer.txt',unpack=True)
    times = np.loadtxt('metalflux_data/'+name[i]+'times.txt',unpack=True)

    dflux = (Zmass_inner - Zmass_outer)/(Zmass_inner + Zmass_outer)
#    print(dflux)
    plt.plot(times,dflux,label=labels[i],color=colors[i])

#plt.text(1,0.6,name,size=15)
plt.text(4,0.7,'More metals within 10 kpc')
plt.text(3.9,-1.11,'More metals outside of 10 kpc')
rect_in  = pat.Rectangle((0,0), 14, 0.8,color='SkyBlue',alpha=0.3)
rect_out = pat.Rectangle((0,-1.201), 14, 1.2,color='Salmon',alpha=0.3)
plt.gca().add_patch(rect_in)
plt.gca().add_patch(rect_out)
plt.ylabel(r'log Metal Flux (at 10 kpc)',fontsize=15)
plt.xlabel('Age/Gyr',fontsize=15)
plt.ylim(-1.2,0.8)
plt.xlim(0,14)
plt.legend()
plt.savefig('dflux_thrutime.pdf')
plt.show()
