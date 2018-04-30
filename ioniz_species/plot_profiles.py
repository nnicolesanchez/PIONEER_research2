import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

labels = ['P0','GM1','GM4','GM5','GM6','GM7']
colors = ['DodgerBlue','SteelBlue','FireBrick','IndianRed','Salmon','Orange']


for k in range(len(labels)-1):
#    T = np.loadtxt(labels[k]+'_T_z0.np')
#    R = np.loadtxt(labels[k]+'_Rbins_z0.np')
    T = np.loadtxt('T_thrutime/'+labels[k]+'_T_4096.np')
    R = np.loadtxt('Rbins_thrutime/'+labels[k]+'_Rbins_4096.np')
    plt.plot(R,T,label=labels[k],color=colors[k])


plt.title('z = 0.0')
plt.ylabel('log(T) [K]')
plt.xlabel('R [kpc]')
plt.ylim(5.2,6.2)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_T_R.pdf')
plt.show()
plt.close()

for j in range(len(labels)-1):
    Novi = np.loadtxt('Novi_thrutime/'+labels[j]+'_Novi_4096.np')
    R = np.loadtxt('Rbins_thrutime/'+labels[j]+'_Rbins_4096.np')
    plt.plot(R,Novi,label=labels[j],color=colors[j])

plt.title('z = 0.0')
plt.ylabel(r'log(N$_{ovi}$) [cm$^{-2}$]')
plt.xlabel('R [kpc]')
plt.ylim(12.5,16.5)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_Novi_R.pdf')
plt.show()
plt.close()

for j in range(len(labels)-1):
    M_ox = np.loadtxt('Omass_thrutime/'+labels[j]+'_Omass_4096.np')
    R = np.loadtxt('Rbins_thrutime/'+labels[j]+'_Rbins_4096.np')
    plt.plot(R,M_ox,label=labels[j],color=colors[j])

plt.title('z = 0.0')
plt.ylabel(r'log(M$_[O]$')
plt.xlabel('R [kpc]')
#plt.ylim(12.5,16.5)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_Omass_R.pdf')
plt.show()
plt.close()

# Calculate cooling time
m_H = 1.67 * 10**-24 #g
C_1 = 3.88 * 10**11 # s K^-1/2 cm^-3
C_2 = 5 * 10**7 #K
f_m = 1.0 # metallicity dependent constant is 1 for solar metallicity Balogh et al 2013
mu = 0.6  # solar metallicity

for j in range(len(labels)-1):
    rho = np.loadtxt('rho_thrutime/'+labels[j]+'_rho_4096.np')
    R = np.loadtxt('Rbins_thrutime/'+labels[j]+'_Rbins_4096.np')
    T = np.loadtxt('T_thrutime/'+labels[k]+'_T_4096.np')
    t_cool = (C_1*mu*m_H*(10**T)**(0.5))/((10**rho)*(1+(C_2*f_m/(10**T)))) / (3.154 * 10**7)# years
    plt.plot(R,np.log10(t_cool),label=labels[j],color=colors[j])

plt.title('z = 0.0')
plt.ylabel(r't$_{cool}$ [years]')
plt.xlabel('R [kpc]')
#plt.ylim(-0.05,1.2)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_tcool_R.pdf')
plt.show()
