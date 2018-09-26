import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

labels = ['P0','GM1','GM7','GM4']
lab_NEW = ['P0','GM1','GM2','GM3'] 
lab_NB = ['P0noBH','GM1noBH','GM7noBH','GM4noBH']
#colors = ['DodgerBlue','SteelBlue','FireBrick','Red','Salmon','Orange']
colors = ['Black','Black','Grey','Grey']
m_p = 1.6726 * 10**-24 #g 

time = '3456'
Z_sun = 0.0142 # (Asplund 2009; https://arxiv.org/pdf/0909.0948.pdf) 
solid = [-0.1,-0.1]
dashed = [-0.1,-0.1]
plt.plot(solid,dashed,color='Black',linestyle='-',label='BH')
plt.plot(solid,dashed,color='Black',linestyle='--',label='NO BH')
plt.plot(solid,dashed,color='Black',linestyle='-',label='Star Forming')
plt.plot(solid,dashed,color='Grey',linestyle='-',label='Quenched')
for k in range(len(labels)):
    Z = np.loadtxt(labels[k]+'_metals_'+time+'.np')
    R = np.loadtxt(labels[k]+'_Rbins_'+time+'.np')
    plt.plot(R,Z/Z_sun,color=colors[k])

for i in range(len(labels)):
    Z_noBH = np.loadtxt(lab_NB[i]+'_metals_3456.np')
    R_noBH = np.loadtxt(lab_NB[i]+'_Rbins_3456.np')
    plt.plot(R_noBH,Z_noBH/Z_sun,linestyle='--',color=colors[0]) #since noBHs are all SF

plt.title('Disk Profile at z = 0.17')
plt.ylabel(r'$Z/Z_{\odot}$',size=15)
plt.xlabel('R [kpc]',size=15)
plt.ylim(0,4)
plt.xlim(0,10)
plt.legend(ncol=2,fontsize=15)
plt.savefig('ALLGMs_DISKmetals_R_plusnoBH.pdf')
plt.show()
plt.close()


quit()









for k in range(len(labels)):
    totgasmass = np.loadtxt('../ioniz_species/totgasmass_thrutime/'+labels[k]+'_totmassgas_'+time+'.np')
    totgasmass = totgasmass/(5.976*10**27) #solar mass
    R = np.loadtxt('../ioniz_species/Rbins_thrutime/'+labels[k]+'_Rbins_'+time+'.np')
    plt.plot(R,np.log10(totgasmass),label=labels[k],color=colors[k])

totgas_noBH = np.loadtxt('GM4noBHs_totgasmass_3456.np')
totgas_noBH = totgas_noBH/(5.976*10**27) #solar mass
R_noBH = np.loadtxt('GM4noBHs_Rbins_3456.np')
plt.plot(R_noBH,np.log10(totgas_noBH),label='GM4_noBH',color='Green',linestyle='--')

plt.title('z = 0.17')
plt.ylabel('log(M$_{gas}$) [M$_{\odot}$]')
plt.xlabel('R [kpc]')
#plt.ylim(5.2,6.2)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_totgasmass_R_plusGM4noBH.pdf')
plt.show()
plt.close()






for k in range(len(labels)):
    T = np.loadtxt('../ioniz_species/T_thrutime/'+labels[k]+'_T_'+time+'.np')
    R = np.loadtxt('../ioniz_species/Rbins_thrutime/'+labels[k]+'_Rbins_'+time+'.np')
    plt.plot(R,T,label=labels[k],color=colors[k])

T_noBH = np.loadtxt('GM4noBHs_T_3456.np')
R_noBH = np.loadtxt('GM4noBHs_Rbins_3456.np')
plt.plot(R_noBH,T_noBH,label='GM4_noBH',color='Green',linestyle='--')

plt.title('z = 0.17')
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
plt.plot(R_noBH,Novi_noBH,label='GM4_noBH',color='Green',linestyle='--')

plt.title('z = 0.17')
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
    M_ox = np.log10((10**M_ox)*16*m_p/(5.97*10**27))
    R = np.loadtxt('../ioniz_species/Rbins_thrutime/'+labels[j]+'_Rbins_'+time+'.np')
    plt.plot(R,M_ox,label=labels[j],color=colors[j])

Mox_noBH = np.loadtxt('GM4noBHs_Omass_3456.np')
Mox_noBH = np.log10((10**Mox_noBH)*16*m_p/(5.97*10**27))
R_noBH = np.loadtxt('GM4noBHs_Rbins_3456.np')
plt.plot(R_noBH,Mox_noBH,label='GM4_noBH',color='Green',linestyle='--')

plt.title('z = 0.17')
plt.ylabel(r'log(M$_{O}$) [M$_{\odot}$]')
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

rho_noBH = np.loadtxt('GM4noBHs_rho_3456.np')
R_noBH = np.loadtxt('GM4noBHs_Rbins_3456.np')
plt.plot(R_noBH,rho_noBH,label='GM4_noBH',color='Green',linestyle='--')

plt.title('z = 0.17')
plt.ylabel(r'log($\rho$) [g cm$^{-3}$]')
plt.xlabel('R [kpc]')
#plt.ylim(12.5,16.5) 
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_rho_R_plusGM4noBH.pdf')
plt.show()
plt.close()


#quit()
# Calculate cooling time
m_H = 1.67 * 10**-24 #g
C_1 = 3.88 * 10**11 # s K^-1/2 cm^-3
C_2 = 5 * 10**7 #K
f_m = 1.0 # metallicity dependent constant is 1 for solar metallicity Balogh et al 2013
mu = 0.6  # solar metallicity

for i in range(len(labels)):
    rho = np.loadtxt('../ioniz_species/rho_thrutime/'+labels[i]+'_rho_'+time+'.np')
    R = np.loadtxt('../ioniz_species/Rbins_thrutime/'+labels[i]+'_Rbins_'+time+'.np')
    T = np.loadtxt('../ioniz_species/T_thrutime/'+labels[i]+'_T_'+time+'.np')
    t_cool = (C_1*mu*m_H*(10**T)**(0.5))/((10**rho)*(1+(C_2*f_m/(10**T)))) / (3.154 * 10**7)# years
    plt.plot(R,np.log10(t_cool),label=labels[i],color=colors[i])

rho_noBH = np.loadtxt('GM4noBHs_rho_'+time+'.np')
R_noBH = np.loadtxt('GM4noBHs_Rbins_'+time+'.np')
T_noBH = np.loadtxt('GM4noBHs_T_'+time+'.np')
t_cool_noBH = (C_1*mu*m_H*(10**T_noBH)**(0.5))/((10**rho_noBH)*(1+(C_2*f_m/(10**T_noBH)))) / (3.154 * 10**7)# years
plt.plot(R,np.log10(t_cool_noBH),label='GM4_noBH',color='Green',linestyle='--')

plt.title('z = 0.17')
plt.ylabel(r'log(t$_{cool}$) [years]')
plt.xlabel('R [kpc]')
#plt.ylim(-0.05,1.2)
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_tcool_R_plusGM4noBH.pdf')
plt.show()

#quit()


for k in range(len(labels)):
    cot = np.loadtxt('../ioniz_species/coolontime_thrutime/'+labels[k]+'_coolontime_'+time+'.np')
    R = np.loadtxt('../ioniz_species/Rbins_thrutime/'+labels[k]+'_Rbins_'+time+'.np')
    plt.plot(R,cot/(3.154*10**7)/(10**9),label=labels[k],color=colors[k])

cot_noBH = np.loadtxt('GM4noBHs_coolontime_3456.np')
R_noBH = np.loadtxt('GM4noBHs_Rbins_3456.np')
plt.plot(R_noBH,cot_noBH/(3.154*10**7)/(10**9),label='GM4_noBH',color='Green',linestyle='--')

plt.title('z = 0.17')
plt.ylabel(r'cool on time [Gyr]')
plt.xlabel('R [kpc]')
#plt.ylim(5.2,6.2)                          
plt.xlim(-10,260)
plt.legend()
plt.savefig('ALLGMs_coolontime_R_plusGM4noBH.pdf')
plt.show()
plt.close()





quit()
