import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt

import numpy as np

labels = ['P0','GM1','GM7','GM4']
NEW_lab = ['P0','GM1','GM2','GM3']
colors = ['DodgerBlue','SteelBlue','FireBrick','Salmon']#''Red','Salmon','Orange']
labels_noBHs = ['P0noBH','GM1noBH','GM7noBH','GM4noBH']
NEW_lab_noBHs = ['P0_noBH','GM1_noBH','GM2_noBH','GM3_noBH']
colors_noBHs = ['DodgerBlue','SteelBlue','FireBrick','Salmon']#,'Orange']
noBHline = '--'

#OLDGM4_noBH = 'OLDGM4noBHs'
#OLDGM4_color = 'Green'
m_p = 1.6726 * 10**-24 #g 
time = '3456'
Z_sun = 0.0142 # (Asplund 2009; https://arxiv.org/pdf/0909.0948.pdf) 

solid = [0,0.1]
dashed = [0,0.1]
plt.plot(solid,dashed,color='Black',linestyle='-',label='BH')
plt.plot(solid,dashed,color='Black',linestyle='--',label='NO BH')
for k in range(len(labels)):
    print('Read in: ',labels[k])
    totgasmass = np.loadtxt(labels[k]+'/'+labels[k]+'_totgasmass_'+time+'.np')
    totgasmass = totgasmass/(5.976*10**27) #solar mass
    R = np.loadtxt(labels[k]+'/'+labels[k]+'_Rbins_'+time+'.np')
    plt.plot(R,np.log10(totgasmass),label=NEW_lab[k],color=colors[k])

for j in range(len(labels_noBHs)):
    print('Read in: ',labels_noBHs[j])
    noBH_totgasmass = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_totgasmass_'+time+'.np')
    noBH_totgasmass = noBH_totgasmass/(5.976*10**27) #solar mass                    
    noBH_R = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_Rbins_'+time+'.np')
    plt.plot(noBH_R,np.log10(noBH_totgasmass),color=colors_noBHs[j],linestyle=noBHline)

plt.title('z = 0.17')
plt.ylabel('log(M$_{gas}$) [M$_{\odot}$]',size=15)
plt.xlabel('R [kpc]',size=15)
plt.ylim(13,15)
plt.xlim(-10,260)
plt.legend(ncol=2,loc=3,fontsize=15)
plt.savefig('ALLGMs_plusnoBH_totgasmass_R.pdf')
plt.show()
plt.close()

solid = [-0.05,-0.1]
dashed = [-0.05,-0.1]
plt.plot(solid,dashed,color='Black',linestyle='-',label='BH')
plt.plot(solid,dashed,color='Black',linestyle='--',label='NO BH')
for k in range(len(labels)):
    Z = np.loadtxt(labels[k]+'/'+labels[k]+'_Z_'+time+'.np')
    R = np.loadtxt(labels[k]+'/'+labels[k]+'_Rbins_'+time+'.np')
    plt.plot(R,Z/Z_sun,label=NEW_lab[k],color=colors[k])

for j in range(len(labels_noBHs)):
    print('Read in: ',labels_noBHs[j])
    noBH_Z = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_Z_'+time+'.np')
    noBH_R = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_Rbins_'+time+'.np')
    plt.plot(noBH_R,noBH_Z/Z_sun,color=colors_noBHs[j],linestyle=noBHline)

#OLDGM4noBH_sfh = np.loadtxt('OLDGM4noBHs_sfh_3456.np')
#OLDGM4noBH_R = np.loadtxt('OLDGM4noBHs_Rbins_3456.np')
#plt.plot(OLDGM4noBH_R,OLDGM4noBH_sfh,label='GM4_noBH',color='Green',linestyle='--')

plt.title('z = 0.17')
plt.ylabel(r'$Z/Z_{\odot}$',size=15)
plt.xlabel('R [kpc]',size=15)
plt.ylim(0,1)
plt.xlim(-10,260)
plt.legend(ncol=2,fontsize=15)
plt.savefig('ALLGMs_plusnoBH_Z_R.pdf')
plt.show()
plt.close()


solid = [0,0.1]
dashed = [0,0.1]
plt.plot(solid,dashed,color='Black',linestyle='-',label='BH')
plt.plot(solid,dashed,color='Black',linestyle='--',label='NO BH')
for k in range(len(labels)):
    T = np.loadtxt(labels[k]+'/'+labels[k]+'_T_'+time+'.np')
    R = np.loadtxt(labels[k]+'/'+labels[k]+'_Rbins_'+time+'.np')
    plt.plot(R,T,label=NEW_lab[k],color=colors[k])

for j in range(len(labels_noBHs)):
    print('Read in: ',labels_noBHs[j])
    noBH_T = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_T_'+time+'.np')
    noBH_R = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_Rbins_'+time+'.np')
    plt.plot(noBH_R,noBH_T,color=colors_noBHs[j],linestyle=noBHline)

#OLDGM4noBH_T = np.loadtxt('OLDGM4noBHs_T_3456.np')
#OLDGM4noBH_R = np.loadtxt('OLDGM4noBHs_Rbins_3456.np')
#plt.plot(OLDGM4noBH_R,OLDGM4noBH_T,label='GM4_noBH',color='Green',linestyle='--')

plt.title('z = 0.17')
plt.ylabel('log(T) [K]',size=15)
plt.xlabel('R [kpc]',size=15)
plt.ylim(5.2,6.2)
plt.xlim(-10,260)
plt.legend(ncol=2,fontsize=15)
plt.savefig('ALLGMs_plusnoBH_T_R.pdf')
plt.show()
plt.close()

solid = [0,0.1]
dashed = [0,0.1]
plt.plot(solid,dashed,color='Black',linestyle='-',label='BH')
plt.plot(solid,dashed,color='Black',linestyle='--',label='NO BH')
for j in range(len(labels)):
    Novi = np.loadtxt(labels[j]+'/'+labels[j]+'_Novi_'+time+'.np')
    R = np.loadtxt(labels[j]+'/'+labels[j]+'_Rbins_'+time+'.np')
    plt.plot(R,Novi,label=NEW_lab[j],color=colors[j])

for j in range(len(labels_noBHs)):
    print('Read in: ',labels_noBHs[j])
    noBH_Novi = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_Novi_'+time+'.np')
    noBH_R = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_Rbins_'+time+'.np')
    plt.plot(noBH_R,noBH_Novi,color=colors_noBHs[j],linestyle=noBHline)

#OLDGM4noBH_Novi = np.loadtxt('OLDGM4noBHs_Novi_3456.np')
#OLDGM4noBH_R = np.loadtxt('OLDGM4noBHs_Rbins_3456.np')
#plt.plot(OLDGM4noBH_R,OLDGM4noBH_Novi,label='GM4_noBH',color='Green',linestyle='--')

plt.title('z = 0.17')
plt.ylabel(r'log(N$_{ovi}$) [cm$^{-2}$]',size=15)
plt.xlabel('R [kpc]',size=15)
plt.ylim(12.5,16.5)
plt.xlim(-10,260)
plt.legend(ncol=2,fontsize=15)
plt.savefig('ALLGMs_plusnoBH_Novi_R.pdf')
plt.show()
plt.close()

solid = [0,0.1]
dashed = [0,0.1]
plt.plot(solid,dashed,color='Black',linestyle='-',label='BH')
plt.plot(solid,dashed,color='Black',linestyle='--',label='NO BH')
for j in range(len(labels)):
    rho = np.loadtxt(labels[j]+'/'+labels[j]+'_rho_'+time+'.np')
    R = np.loadtxt(labels[j]+'/'+labels[j]+'_Rbins_'+time+'.np')
    plt.plot(R,rho,label=NEW_lab[j],color=colors[j])

for j in range(len(labels_noBHs)):
    print('Read in: ',labels_noBHs[j])
    noBH_rho = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_rho_'+time+'.np')
    noBH_R = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_Rbins_'+time+'.np')
    plt.plot(noBH_R,noBH_rho,color=colors_noBHs[j],linestyle=noBHline)

plt.title('z = 0.17')
plt.ylabel(r'log($\rho$) [g cm$^{-3}$]',size=15)
plt.xlabel('R [kpc]',size=15)
plt.ylim(-27,-24.5) 
plt.xlim(-10,260)
plt.legend(fontsize=15,ncol=2)
plt.savefig('ALLGMs_plusnoBH_rho_R.pdf')
plt.show()
plt.close()

