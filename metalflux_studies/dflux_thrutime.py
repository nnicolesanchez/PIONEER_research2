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


def smooth(y, box_pts):
    # moving average smoothing
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


ALL_sims = ['/nobackupp2/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM1.1536gst1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.00','/nobackupp2/nnsanche/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00','/nobackupp2/nnsanche/NO_BHs/pioneer50h243.1536gst1bwK1/pioneer50h243.1536gst1bwK1.00','/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM1.1536gst1bwK1/pioneer50h243GM1.1536gst1bwK1.00','/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM7.1536gst1bwK1/pioneer50h243GM7.1536gst1bwK1.00','/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM4.1536gst1bwK1/pioneer50h243GM4.1536gst1bwK1.00']
name   = ['P0','GM1','GM2','GM3']
labels = ['P0','GM1','GM2','GM3']
colors = ['SteelBlue','DodgerBlue','FireBrick','Salmon']

print('Reading in files from metalflux_data/')
for i in range(len(name)):
    #totmass_inner_t1 = np.loadtxt('metalflux_data/'+name[i]+'totmass_inner.txt',unpack=True)
    #totmass_outer  = np.loadtxt('metalflux_data/'+name[i]+'totmass_outer.txt',unpack=True)
    Zmass_inner = np.loadtxt('metalflux_data/'+name[i]+'Zmass_inner.txt',unpack=True)
    Zmass_outer = np.loadtxt('metalflux_data/'+name[i]+'Zmass_outer.txt',unpack=True)
    times = np.loadtxt('metalflux_data/'+name[i]+'times.txt',unpack=True)

    dt   = []
    dM_z = []
    new_times = []
    for j in range(1,len(times)):
        dt_inyears = (times[j]-times[j-1])*10**9
        dt.append(dt_inyears)
        new_times.append(times[j-1] + ((times[j]-times[j-1])/2))
        dM_t1 = Zmass_inner[j-1] - Zmass_outer[j-1]
        dM_t2 = Zmass_inner[j] - Zmass_outer[j]
        dM_z.append(dM_t2 - dM_t1)


    dM_z_dt = np.array(dM_z)/np.array(dt)
    print(dM_z_dt)
    plt.plot(np.array(new_times),dM_z_dt,label=labels[i],color=colors[i])
#    plt.title('Metal Flow Rate vs Time')
#    plt.show()


import tangos
tangos.all_simulations()
import pylab as p
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

P0_halo = tangos.get_halo("pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.004096/halo_1")
P0_SFR = P0_halo["SFR_histogram"]
P0_BH_accrate = P0_halo.calculate('BH.BH_mdot_histogram')
P0_BH_accrate_hat = smooth(P0_BH_accrate, 50)
P0_SFR_property_object = P0_halo.get_objects("SFR_histogram")[0]
P0_SFR_time_bins = P0_SFR_property_object.x_values()

GM1_halo = tangos.get_halo("pioneer50h243GM1.1536gst1bwK1BH_no3072/pioneer50h243GM1.1536gst1bwK1BH.004096/halo_1")
GM1_SFR = GM1_halo["SFR_histogram"]
GM1_BH_accrate = GM1_halo.calculate('BH.BH_mdot_histogram')
print(GM1_BH_accrate[0:100])
print(np.sum(GM1_BH_accrate))
GM1_BH_accrate_hat = smooth(GM1_BH_accrate, 50)
GM1_SFR_property_object = GM1_halo.get_objects("SFR_histogram")[0]
GM1_SFR_time_bins = GM1_SFR_property_object.x_values()

GM7_halo = tangos.get_halo("pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.003968/halo_1")
GM7_SFR = GM7_halo["SFR_histogram"]
GM7_BH_accrate = GM7_halo.calculate('BH.BH_mdot_histogram')
GM7_BH_accrate_hat = smooth(GM7_BH_accrate, 50)
GM7_SFR_property_object = GM7_halo.get_objects("SFR_histogram")[0]
GM7_SFR_time_bins = GM7_SFR_property_object.x_values()

GM4_halo = tangos.get_halo("pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.004096/halo_1")
GM4_SFR = GM4_halo["SFR_histogram"]
GM4_BH_accrate = GM4_halo.calculate('BH.BH_mdot_histogram')
GM4_BH_accrate_hat = smooth(GM4_BH_accrate, 50)
GM4_SFR_property_object = GM4_halo.get_objects("SFR_histogram")[0]
GM4_SFR_time_bins = GM4_SFR_property_object.x_values()
plt.plot(P0_SFR_time_bins, np.log10(P0_BH_accrate_hat),color='DodgerBlue')#,label='P0')
plt.plot(GM1_SFR_time_bins, np.log10(GM1_BH_accrate_hat),color='SteelBlue')#,label='GM1')
plt.plot(GM7_SFR_time_bins, np.log10(GM7_BH_accrate_hat),color='FireBrick')#,label='GM2')
plt.plot(GM4_SFR_time_bins, np.log10(GM4_BH_accrate_hat),color='Salmon')#,label='GM3')

plt.plot([-10,-11],[-10,-11],color='Black',linestyle='-',label=r'Metal Flow Rate (at 10 kpc) [M$_{\odot}$/yr]')
plt.plot([-10,-11],[-10,-11],color='Black',linestyle='--',label=r'log BH Accretion Rate [M$_{\odot}$ yr$^{-1}$]')


#plt.text(1,0.6,name,size=15)
#plt.text(4,0.7,'More metals within 10 kpc')
#plt.text(3.9,-1.11,'More metals outside of 10 kpc')
#rect_in  = pat.Rectangle((0,0), 14, 0.8,color='SkyBlue',alpha=0.3)
#rect_out = pat.Rectangle((0,-1.201), 14, 1.2,color='Salmon',alpha=0.3)
#plt.gca().add_patch(rect_in)
#plt.gca().add_patch(rect_out)
#plt.ylabel(r'Metal Flow Rate (at 10 kpc) [M$_{\odot}$/yr]',fontsize=15)
plt.ylabel('Rates [M$_{\odot}$/yr]')
plt.xlabel('Age/Gyr',fontsize=15)
#plt.yscale('log')
plt.ylim(-6,1)
plt.xlim(0,14)
plt.legend(ncol=2)
plt.savefig('metalflowrate_plusBHacrr_thrutime.pdf')
plt.show()

