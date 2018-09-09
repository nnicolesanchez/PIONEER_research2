import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt

import numpy as np

labels = ['P0','GM1','GM7','GM4']
NEW_lab = ['P0','GM1','GM2','GM3']
colors = ['DodgerBlue','SteelBlue','FireBrick','Salmon']#''Red','Salmon','Orange']
labels_noBHs = ['P0noBH','GM1noBH','GM7noBH','GM4noBH']
NEW_lab_noBHs = ['P0_noBH','GM1_noBH','GM7_noBH','GM4_noBH']
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
    totgasmass = totgasmass # Converted arrays to Msun /(5.976*10**27) #solar mass
    R = np.loadtxt(labels[k]+'/'+labels[k]+'_Rbins_'+time+'.np')
    plt.plot(R,np.log10(totgasmass),label=NEW_lab[k],color=colors[k])

for j in range(len(labels_noBHs)):
    print('Read in: ',labels_noBHs[j])
    noBH_totgasmass = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_totgasmass_'+time+'.np')
    noBH_totgasmass = noBH_totgasmass# Converted arrays to Msun /(5.976*10**27) #solar mass                    
    noBH_R = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_Rbins_'+time+'.np')
    plt.plot(noBH_R,np.log10(noBH_totgasmass),color=colors_noBHs[j],linestyle=noBHline)

plt.title('z = 0.17')
plt.ylabel('log(M$_{gas}$) [M$_{\odot}$]',size=15)
plt.xlabel('R [kpc]',size=15)
plt.ylim(7.6,9.6)
plt.xlim(-10,260)
#plt.legend(ncol=2,loc=3,fontsize=15)
plt.savefig('ALLGMs_plusnoBH_totgasmass_R.pdf')
plt.show()
plt.close()


solid = [0,0.1]
dashed = [0,0.1]
plt.plot(solid,dashed,color='Black',linestyle='-',label='BH')
plt.plot(solid,dashed,color='Black',linestyle='--',label='NO BH')
for k in range(len(labels)):
    print('Read in: ',labels[k])
    Omass = np.loadtxt(labels[k]+'/'+labels[k]+'_Omass_'+time+'.np')
    Omass = Omass # Converted arrays to Msun /(5.976*10**27) #solar mass                                                  
    R = np.loadtxt(labels[k]+'/'+labels[k]+'_Rbins_'+time+'.np')
    plt.plot(R,np.log10(Omass),label=NEW_lab[k],color=colors[k])

for j in range(len(labels_noBHs)):
    print('Read in: ',labels_noBHs[j])
    noBH_Omass = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_Omass_'+time+'.np')
    noBH_Omass = noBH_Omass# Converted arrays to Msun /(5.976*10**27) #solar mass                                         
    noBH_R = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_Rbins_'+time+'.np')
    plt.plot(noBH_R,np.log10(noBH_Omass),color=colors_noBHs[j],linestyle=noBHline)

plt.title('z = 0.17')
plt.ylabel('log(M$_{oxygen gas}$) [M$_{\odot}$]',size=15)
plt.xlabel('R [kpc]',size=15)
#plt.ylim(7.6,9.6)
plt.xlim(-10,260)
#plt.legend(ncol=2,loc=3,fontsize=15)                                                                                              
plt.savefig('ALLGMs_plusnoBH_Omass_R.pdf')
plt.show()
plt.close()


#solid = [-0.05,-0.1]
#dashed = [-0.05,-0.1]
#plt.plot(solid,dashed,color='Black',linestyle='-',label='BH')
#plt.plot(solid,dashed,color='Black',linestyle='--',label='NO BH')
for k in range(len(labels)):
    Z = np.loadtxt(labels[k]+'/'+labels[k]+'_Z_'+time+'.np')
    R = np.loadtxt(labels[k]+'/'+labels[k]+'_Rbins_'+time+'.np')
    plt.plot(R,np.log10(Z/Z_sun),label=NEW_lab[k],color=colors[k])

for j in range(len(labels_noBHs)):
    print('Read in: ',labels_noBHs[j])
    noBH_Z = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_Z_'+time+'.np')
    noBH_R = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_Rbins_'+time+'.np')
    plt.plot(noBH_R,np.log10(noBH_Z/Z_sun),color=colors_noBHs[j],linestyle=noBHline)

#OLDGM4noBH_sfh = np.loadtxt('OLDGM4noBHs_sfh_3456.np')
#OLDGM4noBH_R = np.loadtxt('OLDGM4noBHs_Rbins_3456.np')
#plt.plot(OLDGM4noBH_R,OLDGM4noBH_sfh,label='GM4_noBH',color='Green',linestyle='--')

plt.title('z = 0.17')
plt.ylabel(r'log($Z/Z_{\odot}$)',size=15)
plt.xlabel('R [kpc]',size=15)
plt.ylim(-1.6,0.5)
plt.xlim(-10,260)
#plt.legend(ncol=2,fontsize=15)
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
#plt.show()
plt.close()

solid = [0,0.1]
dashed = [0,0.1]
plt.plot(solid,dashed,color='Black',linestyle='-',label='BH')
plt.plot(solid,dashed,color='Black',linestyle='--',label='NO BH')
for j in range(len(labels)):
    Novi = np.loadtxt(labels[k]+'/'+labels[k]+'_Novi_'+time+'.np')
    R = np.loadtxt(labels[k]+'/'+labels[k]+'_Rbins_'+time+'.np')
#    Novi = np.loadtxt('/home1/nnsanche/PIONEER_research2/ioniz_species/'+labels[j]+'_Novi_'+time+'_hdf5.np')
#    R = np.loadtxt('/home1/nnsanche/PIONEER_research2/ioniz_species/'+labels[j]+'_Rbins_'+time+'.np')
    plt.plot(R,Novi,label=NEW_lab[j],color=colors[j])

for j in range(len(labels_noBHs)):
    print('Read in: ',labels_noBHs[j])
    noBH_Novi = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_Novi_'+time+'.np')
    noBH_R = np.loadtxt(labels_noBHs[j]+'/'+labels_noBHs[j]+'_Rbins_'+time+'.np')
#    noBH_Novi = np.loadtxt('/home1/nnsanche/PIONEER_research2/noBH_analysis/'+NEW_lab_noBHs[j]+'_Novi_'+time+'_hdf5.np')
#    noBH_R = np.loadtxt('/home1/nnsanche/PIONEER_research2/noBH_analysis/'+NEW_lab_noBHs[j]+'_Rbins_'+time+'.np')
    plt.plot(noBH_R,noBH_Novi,color=colors_noBHs[j],linestyle=noBHline)

COS = pd.read_csv('COSHALOdata.txt',comment='#',index_col=False,header=None)

print('Average z for COS halos:',COS[1].mean())
print(COS[2],COS[15])
for p in range(len(COS[0])):
    print(COS[0][p])
    if COS[17][p] == ' RED' :
        mark = 's'
        markerface = 'IndianRed'
        lab = 'COS Ellipticals'
    else:
        mark = 'o'
        markerface = 'SteelBlue'
        lab = 'COS Spirals'
    if COS[14][p] == ' <':
        upperlims = True
        lowerlims = False
        errors = 0.1
        markerface = 'none'
    elif COS[14][p] == ' >':
        upperlims = False
        lowerlims = True
        errors = 0.1
    else:
        upperlims = False
        lowerlims = False
        errors = 0

    print(COS[2][p],COS[15][p])
    plt.errorbar(COS[2][p],COS[15][p],yerr=errors,markerfacecolor=markerface,marker=mark,markersize=7,
color='Black',uplims=upperlims,lolims=lowerlims)

plt.plot([1,2],[1,2],marker='s',label='COS Elliptical',color='Black',mfc='IndianRed',linestyle='None')
plt.plot([1,2],[1,2],marker='o',label='COS Spiral',color='Black',mfc='SteelBlue',linestyle='None')


plt.title('z = 0.17')
plt.ylabel(r'log(N$_{ovi}$) [cm$^{-2}$]',size=15)
plt.xlabel('R [kpc]',size=14)
plt.ylim(12.5,16.5)
plt.xlim(-10,260)
plt.legend(ncol=2,fontsize=15)
plt.savefig('ALLGMs_plusnoBH_Novi_R.pdf')
plt.show()
plt.close()

#### NOVII PROFILE ####
solid = [0,0.1]
dashed = [0,0.1]
plt.plot(solid,dashed,color='Black',linestyle='-',label='BH')
plt.plot(solid,dashed,color='Black',linestyle='--',label='NO BH')
for j in range(len(labels)):
    Novii = np.loadtxt('/home1/nnsanche/PIONEER_research2/ioniz_species/'+labels[j]+'_Novii_'+time+'_hdf5.np')
    R = np.loadtxt('/home1/nnsanche/PIONEER_research2/ioniz_species/'+labels[j]+'_Rbins_'+time+'.np')
    plt.plot(R,Novii,label=NEW_lab[j],color=colors[j])

for j in range(len(labels_noBHs)):
    print('Read in: ',labels_noBHs[j])
    noBH_Novii = np.loadtxt('/home1/nnsanche/PIONEER_research2/noBH_analysis/'+NEW_lab_noBHs[j]+'_Novii_'+time+'_hdf5.np')
    noBH_R = np.loadtxt('/home1/nnsanche/PIONEER_research2/noBH_analysis/'+NEW_lab_noBHs[j]+'_Rbins_'+time+'.np')
    plt.plot(noBH_R,noBH_Novii,color=colors_noBHs[j],linestyle=noBHline)
plt.title('z = 0.17')
plt.ylabel(r'log(N$_{ovii}$) [cm$^{-2}$]',size=15)
plt.xlabel('R [kpc]',size=14)
#plt.ylim(12.5,16.5)
plt.xlim(-10,260)
plt.legend(ncol=2,fontsize=15)
plt.savefig('ALLGMs_plusnoBH_Novii_R.pdf')
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
plt.ylim(-28.7,-24.3) 
plt.xlim(-10,260)
#plt.legend(fontsize=15,ncol=2)
plt.savefig('ALLGMs_plusnoBH_rho_R.pdf')
plt.show()
plt.close()

