# This script uses the tangos database information to build
# a merger tree for the main halos in the GM sims
import tangos
from tangos.examples import mergers
import numpy as np
import pylab as p

# Main Halos at z = 0
P0_halo1_z0 = tangos.get_halo("pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.004096/halo_1")
print(mergers.get_mergers_of_major_progenitor(P0_halo1_z0))  # Output: redshift, ratio, halos
GM1_halo1_z0 = tangos.get_halo("pioneer50h243GM1.1536gst1bwK1BH_no3072/pioneer50h243GM1.1536gst1bwK1BH.004096/halo_1")
print(mergers.get_mergers_of_major_progenitor(GM1_halo1_z0))
GM4_halo1_z0 = tangos.get_halo("pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.004096/halo_1")
print(mergers.get_mergers_of_major_progenitor(GM4_halo1_z0))
GM7_halo1_z0 = tangos.get_halo("pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.003968/halo_1")
print(mergers.get_mergers_of_major_progenitor(GM7_halo1_z0))

# Plot Halo Mass as function of Time
P0_halo1_Mvir, P0_halo1_t = P0_halo1_z0.calculate_for_progenitors("Mvir","t()")
GM1_halo1_Mvir, GM1_halo1_t = GM1_halo1_z0.calculate_for_progenitors("Mvir","t()")
GM7_halo1_Mvir, GM7_halo1_t = GM7_halo1_z0.calculate_for_progenitors("Mvir","t()")
GM4_halo1_Mvir, GM4_halo1_t = GM4_halo1_z0.calculate_for_progenitors("Mvir","t()")
p.plot(P0_halo1_t,np.log10(P0_halo1_Mvir),color='DodgerBlue',label='P0')
p.plot(GM1_halo1_t,np.log10(GM1_halo1_Mvir),color='SteelBlue',label='GM1')
p.plot(GM7_halo1_t,np.log10(GM7_halo1_Mvir),color='FireBrick',label='GM2')
p.plot(GM4_halo1_t,np.log10(GM4_halo1_Mvir),color='Salmon',label='GM3')
p.xlabel("t/Gyr")
p.ylabel(r"$M/M_{\odot}$")
p.legend()
p.savefig('ALLGM_MHgrowth.pdf')
p.show()


# Satellite Halos at z ~ 1
GM4_halo1_z0 = tangos.get_halo("pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.004096/halo_1")
mergers.get_mergers_of_major_progenitor(GM4_halo1_z0)
GM4_halo3_z1 = tangos.get_halo("pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.001792/halo_3")
GM4_halo1_z1 = tangos.get_halo("pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.001792/halo_1")
GM4_halo3_z1['Mvir']
GM4_halo1_z1['Mvir']





