import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.stats import *

Mgas     = np.loadtxt('COSstellarmass_range_data/ROM_Mtotgas_10_100kpc_z017.txt')
Mmetals  = np.loadtxt('COSstellarmass_range_data/ROM_Mtotmetals_10_100kpc_z017.txt')
Mhalo    = np.loadtxt('COSstellarmass_range_data/ROM_halomasses_z017.txt')
Mstellar = np.loadtxt('COSstellarmass_range_data/ROM_stellarmass_z017.txt')

BHmass = np.loadtxt('ROM_data_forhalos/ROM_BH_masses.txt')

Mbaryon         = Mhalo*0.17
Mgas_Mbaryon    = Mgas/Mbaryon
Mmetals_Mbaryon = Mmetals/Mbaryon

rho, p_val = spearmanr(BHmass,Mgas)
print(rho,p_val)

# Plot Mgas/Mbaryon vs BHmass
sm = plt.scatter(np.log10(BHmass),Mgas_Mbaryon,marker='o',alpha=0.7,cmap='rainbow',c=np.log10(Mstellar),vmin=10.5,vmax=11.5)

plt.ylabel(r'M$_{gas}$/M$_{baryon}$',size=15)
plt.xlabel(r'M$_{BH}$/M$_{\odot}$',size=15)
plt.colorbar().set_label(label=r'log M$_{*}$/M$_{\odot}$',size=15)
plt.savefig('Mgasoverbaryons_BHmass.pdf')
plt.show()

