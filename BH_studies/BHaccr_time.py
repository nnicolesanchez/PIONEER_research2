#import matplotlib
#matplotlib.use('Agg')

import tangos
tangos.all_simulations()
import pylab as p
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

P0_halo = tangos.get_halo("pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.004096/halo_1")
P0_SFR = P0_halo["SFR_histogram"]
P0_BH_accrate = P0_halo.calculate('BH.BH_mdot_histogram')
P0_SFR_property_object = P0_halo.get_objects("SFR_histogram")[0]
P0_SFR_time_bins = P0_SFR_property_object.x_values()

GM1_halo = tangos.get_halo("pioneer50h243GM1.1536gst1bwK1BH_no3072/pioneer50h243GM1.1536gst1bwK1BH.004096/halo_1")
GM1_SFR = GM1_halo["SFR_histogram"]
GM1_BH_accrate = GM1_halo.calculate('BH.BH_mdot_histogram')
GM1_SFR_property_object = GM1_halo.get_objects("SFR_histogram")[0]
GM1_SFR_time_bins = GM1_SFR_property_object.x_values()

GM7_halo = tangos.get_halo("pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.003968/halo_1")
GM7_SFR = GM7_halo["SFR_histogram"]
GM7_BH_accrate = GM7_halo.calculate('BH.BH_mdot_histogram')
GM7_SFR_property_object = GM7_halo.get_objects("SFR_histogram")[0]
GM7_SFR_time_bins = GM7_SFR_property_object.x_values()

GM4_halo = tangos.get_halo("pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.004096/halo_1")
GM4_SFR = GM4_halo["SFR_histogram"]
GM4_BH_accrate = GM4_halo.calculate('BH.BH_mdot_histogram')
GM4_SFR_property_object = GM4_halo.get_objects("SFR_histogram")[0]
GM4_SFR_time_bins = GM4_SFR_property_object.x_values()


plt.plot(P0_SFR_time_bins, np.log10(P0_BH_accrate),color='DodgerBlue',label='P0')
plt.plot(GM1_SFR_time_bins, np.log10(GM1_BH_accrate),color='SteelBlue',label='GM1')
plt.plot(GM7_SFR_time_bins, np.log10(GM7_BH_accrate),color='FireBrick',label='GM2')
plt.plot(GM4_SFR_time_bins, np.log10(GM4_BH_accrate),color='Salmon',label='GM3')
plt.xlabel("Age/Gyr")
plt.ylabel("BH accretion rate/$M_{\odot}\,yr^{-1}$")
plt.ylim(-6,-1)
plt.legend()
plt.savefig('ALL_bhaccr_age.pdf')
#plt.show()
plt.clf()

P0_BH = P0_halo['BH'][0]
GM1_BH = GM1_halo['BH'][0]
GM7_BH = GM7_halo['BH'][0]
GM4_BH = GM4_halo['BH'][0]

print(GM7_BH.keys())
plt.plot(P0_BH.calculate_for_progenitors('t()')[0][:-1],np.log10(P0_BH.calculate_for_progenitors('BH_mass')[0]),color='DodgerBlue',label='P0')
plt.plot(GM1_BH.calculate_for_progenitors('t()')[0][:-3],np.log10(GM1_BH.calculate_for_progenitors('BH_mass')[0]),color='SteelBlue',label='GM1')
plt.plot(GM7_BH.calculate_for_progenitors('t()')[0],np.log10(GM7_BH.calculate_for_progenitors('BH_mass')[0]),color='FireBrick',label='GM2')
plt.plot(GM4_BH.calculate_for_progenitors('t()')[0][:-2],np.log10(GM4_BH.calculate_for_progenitors('BH_mass')[0]),color='Salmon',label='GM3')
smooth_data = pd.rolling_mean(ts,5).plot(style='k')
plt.ylabel(r"log(M$_{BH}$/M$_{\odot}$)")
plt.xlabel("Age/Gyr")
plt.legend()
plt.savefig('ALL_bhmass_age.pdf')
plt.show()
