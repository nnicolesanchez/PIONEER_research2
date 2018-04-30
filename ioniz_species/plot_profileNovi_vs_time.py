#
#
#
#

#
#
import matplotlib.pyplot as plt
import numpy as np

k=0
labels = ['P0','GM1','GM4','GM5','GM6']

#for k in range(len(labels)):
N_ovi_50kpc_thru_time = np.loadtxt(labels[k]+'_Novi50kpcvstime.np')
times = np.loadtxt(labels[k]+'_times.np')

#print(N_ovi_50kpc_thru_time,times)

plt.plot(times,np.log10(N_ovi_50kpc_thru_time),label=labels[k])

plt.ylabel(r'N$_{OVI}$ (cm$^{-2}$)')
plt.xlabel('time [Gyr]')
plt.ylim(12,15.5)
plt.xlim(3,14)
plt.legend()
#plt.savefig('ALL_Novi_profile_vs_time_AHF.pdf')
plt.savefig(labels[k]+'_Novi_profile_vs_time_AHF.pdf')
plt.show()
