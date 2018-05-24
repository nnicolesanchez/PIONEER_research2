import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

labels = ['P0','GM1','GM4','GM5','GM6','GM7']
colors = ['DodgerBlue','SteelBlue','FireBrick','IndianRed','Salmon','Orange']

time = '3456'
for k in range(len(labels)):
    T = np.loadtxt('../ioniz_species/T_thrutime/'+labels[k]+'_T_'+time+'.np')
    R = np.loadtxt('../ioniz_species/Rbins_thrutime/'+labels[k]+'_Rbins_'+time+'.np')
    plt.plot(R,T,label=labels[k],color=colors[k])

T_noBH = np.loadtxt('GM4noBHs_T_3456.np')
R_noBH = np.loadtxt('GM4noBHs_Rbins_3456.np')
plt.plot(R_noBH,T_noBH,label='GM4_noBH',color='Green')

plt.title('z = 0.0')
plt.ylabel('log(T) [K]')
plt.xlabel('R [kpc]')
plt.ylim(5.2,6.2)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_T_R_plusGM4noBH.pdf')
plt.show()
plt.close()


for j in range(len(labels)):
    Novi = np.loadtxt('../ioniz_species/Novi_thrutime/'+labels[j]+'_Novi_'+time+'.np')
    R = np.loadtxt('../ioniz_species/Rbins_thrutime/'+labels[j]+'_Rbins_'+time+'.np')
    plt.plot(R,Novi,label=labels[j],color=colors[j])

Novi_noBH = np.loadtxt('GM4noBHs_Novi_3456.np')
R_noBH = np.loadtxt('GM4noBHs_Rbins_3456.np')
plt.plot(R_noBH,Novi_noBH,label='GM4_noBH',color='Green')

plt.title('z = 0.0')
plt.ylabel(r'log(N$_{ovi}$) [cm$^{-2}$]')
plt.xlabel('R [kpc]')
plt.ylim(12.5,16.5)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_Novi_R_plusGM4noBH.pdf')
plt.show()
plt.close()

for j in range(len(labels)):
    M_ox = np.loadtxt('../ioniz_species/Omass_thrutime/'+labels[j]+'_Omass_'+time+'.np')
    R = np.loadtxt('../ioniz_species/Rbins_thrutime/'+labels[j]+'_Rbins_'+time+'.np')
    plt.plot(R,M_ox,label=labels[j],color=colors[j])

Mox_noBH = np.loadtxt('GM4noBHs_Omass_3456.np')
R_noBH = np.loadtxt('GM4noBHs_Rbins_3456.np')
plt.plot(R_noBH,Mox_noBH,label='GM4_noBH',color='Green')

plt.title('z = 0.0')
plt.ylabel(r'log(M$_[O]$')
plt.xlabel('R [kpc]')
#plt.ylim(12.5,16.5)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_Omass_R_plusGM4noBH.pdf')
plt.show()
plt.close()

for j in range(len(labels)):
    rho = np.loadtxt('../ioniz_species/rho_thrutime/'+labels[j]+'_rho_'+time+'.np')
    R = np.loadtxt('../ioniz_species/Rbins_thrutime/'+labels[j]+'_Rbins_'+time+'.np')
    plt.plot(R,rho,label=labels[j],color=colors[j])

Mox_noBH = np.loadtxt('GM4noBHs_rho_3456.np')
R_noBH = np.loadtxt('GM4noBHs_Rbins_3456.np')
plt.plot(R_noBH,Mox_noBH,label='GM4_noBH',color='Green')

plt.title('z = 0.0')
plt.ylabel(r'log($\rho$)')
plt.xlabel('R [kpc]')
#plt.ylim(12.5,16.5) 
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_rho_R_plusGM4noBH.pdf')
plt.show()
plt.close()







quit()
# Calculate cooling time
m_H = 1.67 * 10**-24 #g
C_1 = 3.88 * 10**11 # s K^-1/2 cm^-3
C_2 = 5 * 10**7 #K
f_m = 1.0 # metallicity dependent constant is 1 for solar metallicity Balogh et al 2013
mu = 0.6  # solar metallicity

for i in range(len(labels)):
    if labels[i] == 'GM7':
        time = '3968'
    else:
        time = '4096'

    rho = np.loadtxt('rho_thrutime/'+labels[i]+'_rho_'+time+'.np')
    R = np.loadtxt('Rbins_thrutime/'+labels[i]+'_Rbins_'+time+'.np')
    T = np.loadtxt('T_thrutime/'+labels[i]+'_T_'+time+'.np')
    t_cool = (C_1*mu*m_H*(10**T)**(0.5))/((10**rho)*(1+(C_2*f_m/(10**T)))) / (3.154 * 10**7)# years
    plt.plot(R,np.log10(t_cool),label=labels[i],color=colors[i])

plt.title('z = 0.0')
plt.ylabel(r'log(t$_{cool}$) [years]')
plt.xlabel('R [kpc]')
#plt.ylim(-0.05,1.2)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_tcool_R.pdf')
plt.show()


for k in range(len(labels)):
    if labels[k] == 'GM7':
        time = '3968'
    else:
        time = '4096'

    cot = np.loadtxt('coolontime_thrutime/'+labels[k]+'_coolontime_'+time+'.np')
    R = np.loadtxt('Rbins_thrutime/'+labels[k]+'_Rbins_'+time+'.np')
    plt.plot(R,cot/(3.154*10**7)/(10**9),label=labels[k],color=colors[k])


plt.title('z = 0.0')
plt.ylabel(r'cool on time [Gyr]')
plt.xlabel('R [kpc]')
#plt.ylim(5.2,6.2)                                                                                                                                                                                                 
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_coolontime_R.pdf')
plt.show()
plt.close()
quit()

for k in range(len(labels)):
    if labels[k] == 'GM7':
        time = '3968'
    else:
        time = '4096'

    Z = np.loadtxt('metals_thrutime/'+labels[k]+'_metals_'+time+'.np')
    R = np.loadtxt('Rbins_thrutime/'+labels[k]+'_Rbins_'+time+'.np')
    plt.plot(R,np.log(Z),label=labels[k],color=colors[k])


plt.title('z = 0.0')
plt.ylabel(r'log($Z$)')
plt.xlabel('R [kpc]')
#plt.ylim(5.2,6.2)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_metals_R.pdf')
plt.show()
plt.close()
quit()
