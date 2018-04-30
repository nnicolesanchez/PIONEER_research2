import matplotlib.pyplot as plt
import numpy as np
sfh4,bins4 = np.loadtxt('GM4_sfhistory_bins.txt', unpack=True)
sfh1,bins1 = np.loadtxt('GM1_sfhistory_bins.txt', unpack=True)
plt.plot(bins4,sfh4,drawstyle='steps',color='SteelBlue')
plt.plot(bins1,sfh1,drawstyle='steps',color='Salmon')
plt.xlim(0,14)
plt.ylim(0,12)
plt.xlabel('Time [Gyr]', fontsize='large')
plt.ylabel('SFR [M$_\odot$ yr$^{-1}$]', fontsize='large')
plt.show()

