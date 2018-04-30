import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

labels = ['P0','GM1','GM4','GM5','GM6','GM7']
colors = ['DodgerBlue','SteelBlue','FireBrick','IndianRed','Salmon','Orange']


for k in range(len(labels)):
    T = np.loadtxt(labels[k]+'_T_z0.np')
    R = np.loadtxt(labels[k]+'_Rbins_z0.np')
    plt.plot(R,T,label=labels[k],color=colors[k])


plt.title('z = 0.0')
plt.ylabel('log(T) [K]')
plt.xlabel('R [kpc]')
plt.ylim(5.2,6.2)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_T_R.pdf')
plt.show()


for j in range(len(labels)):
    Novi = np.loadtxt(labels[j]+'_N_OVI_z0.np')
    R = np.loadtxt(labels[j]+'_Rbins_z0.np')
    plt.plot(R,Novi,label=labels[j],color=colors[j])

plt.title('z = 0.0')
plt.ylabel(r'log(N$_{ovi}$) [cm$^{-2}$]')
plt.xlabel('R [kpc]')
plt.ylim(12.5,16.5)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_Novi_R.pdf')
plt.show()
