import matplotlib.pyplot as plt
import numpy as np

Mgas = np.loadtxt('COSstellarmass_range_data/ROM_Mtotgas_10_100kpc_z017.txt')
Mmetals = np.loadtxt('COSstellarmass_range_data/ROM_Mtotmetals_10_100kpc_z017.txt')
Mhalo = np.loadtxt('COSstellarmass_range_data/ROM_halomasses_z017.txt')

BHmass = np.loadtxt('ROM_data_forhalos/ROM_BH_masses.txt')

Mbaryon         = Mhalo*0.17
Mgas_Mbaryon    = Mgas/Mbaryon
Mmetals_Mbaryon = Mmetals/Mbaryon

# Plot Mgas/Mbaryon vs BHmass
plt.plot(np.log10(BHmass),Mgas_Mbaryon,marker='o',linestyle=' ',color='Black',alpha=0.7)
plt.ylabel(r'M$_{gas}$/M$_{baryon}$')
plt.xlabel(r'M$_{BH}$/M$_{\odot}$')
plt.savefig('Mgasoverbaryons_BHmass.pdf')
plt.show()

plt.plot(np.log10(BHmass),Mmetals_Mbaryon,marker='o',linestyle=' ',color='Black',alpha=0.7)
plt.ylabel(r'M$_{metals}$/M$_{baryon}$')
plt.xlabel(r'M$_{BH}$/M$_{\odot}$')
plt.savefig('Mmetalsoverbaryons_BHmass.pdf')
plt.show()

