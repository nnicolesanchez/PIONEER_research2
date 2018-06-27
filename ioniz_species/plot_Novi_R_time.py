import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


k = 2
labels = ['P0','GM1','GM4','GM5','GM6','GM7']
ts = np.loadtxt('../'+labels[k]+'/timesteps.txt',dtype=str)
z = np.loadtxt('../'+labels[k]+'/steps_z_time.txt',delimiter=' ',dtype=str,usecols=2,unpack=True,skiprows=1)
print(labels[k])


for t in range(len(ts)):
    t = len(ts)-1
    print(ts[t],z[t])
    Novi = np.loadtxt(labels[k]+'_ionizNovi_'+ts[t]+'_ydata.np')
    R = np.loadtxt(labels[k]+'_ionizb_'+ts[t]+'_xdata.np')

    if (float(ts[t]) >= 2816.):
        COS = pd.read_csv('COShalo_obs.txt',header=0,delim_whitespace=True,comment='#',index_col=False)
        COS_ID = COS['ID']
        COS_OVI = COS['LogN0VI']
        COS_Rkpc = COS['Rho(kpc)'][COS_OVI != 0]
        COS_OVI = COS_OVI[COS_OVI != 0]
#        print(COS)
        plt.plot(COS_Rkpc[COS['RED?'] == 'yes'],COS_OVI[COS['RED?'] == 'yes'],marker='.',color='Red',linestyle=' ',label='COS Ellipticals')
        plt.plot(COS_Rkpc[COS['RED?'] == 'no'],COS_OVI[COS['RED?'] == 'no'],marker='.',color='DodgerBlue',linestyle=' ',label='COS Spirals')

    plt.plot(R,np.log10(Novi),marker='.',label=labels[k])
    plt.ylabel(r'N$_{OVI}$ (photoionized) [cm$^{-2}$]')
    plt.xlabel(r'$r$ [kpc]')
    plt.ylim(13,16)
    plt.xlim(-10,270)
    plt.title(z[t])
    #plt.text(150,15.75,z[t])
    plt.legend()
    plt.savefig(labels[k]+'_photoNOVI_b_'+ts[t]+'.pdf')
    #plt.show()
    plt.clf()
