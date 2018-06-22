# This script reads in the Tangos Database specified by the .bashrc path
# It creates a:
#       - SFR plot
#       - BH mass plot
#       - BH accretion plot

# N. Nicole Sanchez   -- /nobackupp2/nnsanche/PIONEER_research2/BH_studies/bh_accr_plots.py
# Univ. of W, Seattle -- Created: June 5, 2018
import matplotlib.pyplot as plt
import pylab as p
import tangos

# See Tangos/pynbody documentation for some details about this method of plot making with Tangos
P0_halo = tangos.get_halo("pioneer50h243.1536g1bwK1BH/%4096/halo_1")
#GM4_halo = tangos.get_halo("pioneer50h243GM4.1536gst1bwK1BH/%4096/halo_1")
GM7_halo = tangos.get_halo("pioneer50h243GM7.1536gst1bwK1BH/%3968/halo_1")
print(P0_halo.keys(),GM7_halo.keys())

print('Currently reading in simulation and halo',P0_halo)
print('Currently reading in simulation and halo',GM7_halo)

print('First: Sanity check.')
print('Does the SFR plot look like it should compared to SFH?')
P0_SFR = P0_halo["SFR_histogram"]
GM7_SFR = GM7_halo["SFR_histogram"]
# From Tangos DOC: "The above is sufficient to retrieve the histogram; however you probably also want to check
# the size of the time bins. The easiest approach is to request a suitable time array to go with
# the SF history:"
P0_SFR_property_object = P0_halo.get_objects("SFR_histogram")[0]
P0_SFR_time_bins = P0_SFR_property_object.x_values()
p.plot(P0_SFR_time_bins, P0_SFR)

GM7_SFR_property_object = GM7_halo.get_objects("SFR_histogram")[0]
GM7_SFR_time_bins = GM7_SFR_property_object.x_values()
p.plot(GM7_SFR_time_bins, GM7_SFR)

p.xlabel("Time/Gyr")
p.ylabel("SFR/$M_{\odot}\,yr^{-1}$")
plt.savefig('ALL_SFR.pdf')
p.show()

print('Next: Make a plot of the BH mass as a function of time.')
P0_BH_mass_range, P0_t, P0_z = P0_halo.calculate_for_progenitors("BH_central.BH_mass","t()","z()")
plt.plot(P0_t,P0_BH_mass_range)

#GM7_BH_mass_range, GM7_t, GM7_z = GM7_halo.calculate_for_progenitors("BH_central.BH_mass","t()","z()")
#plt.plot(GM7_t,GM7_BH_mass_range)
p.xlabel("Time/Gyr")
p.ylabel("M$_{BH}$/$M_{\odot}$")
plt.savefig('ALL_BHmass.pdf')
p.show()

P0_BH_mass = P0_halo.calculate('BH_central.BH_mass')
P0_BH_max_mass = P0_halo.calculate('link(BH_central, BH_mass, "max")')
P0_BH_closest = P0_halo.calculate('link(BH_central, BH_central_distance, "min")')

#GM7_BH_mass = GM7_halo.calculate('BH_central.BH_mass')
#GM7_BH_max_mass = GM7_halo.calculate('link(BH_central, BH_mass, "max")')
#GM7_BH_closest = GM7_halo.calculate('link(BH_central, BH_central_distance, "min")')
#print('Is the BH with the max mass also the closest to the center?')
#P0_BH_max_mass == P0_BH_closes
#GM7_BH_max_mass == GM7_BH_closest

print('Finally: Plot BH accretion rate as a function of time.')
P0_BH_accrate = P0_halo.calculate('BH.BH_mdot_histogram')
p.plot(P0_SFR_time_bins, P0_BH_accrate)

GM7_BH_accrate = GM7_halo.calculate('BH.BH_mdot_histogram')
p.plot(GM7_SFR_time_bins, GM7_BH_accrate)
p.xlabel("Time/Gyr")
p.ylabel("BH accretion rate/$M_{\odot}\,yr^{-1}$")
plt.savefig('ALL_BHaccr_rate.pdf')
p.show()


