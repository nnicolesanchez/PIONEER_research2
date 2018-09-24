#import matplotlib
#matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from hdf5_script import *
import sys

#catalogue = 'grp'
catalogue = 'AHF'

HALO_data = pd.read_csv('/home1/nnsanche/ROMULUS_research/OVI_columndensities/MWmass_data/ROMULUS_highmass_info_ALL.csv',header=0)

if sys.argv[1] == 'COS':
    #### COS-HALOS Range
    print('Plotting COS-Halos stellas mass range of ROM galaxies')    

    M_blue = [9,10,14,15,16,17,18,19,20,21,22,26,27,28,29,30,31,32,33,34,35,36,37,38,42,43,44,45,46,47,48,49,50,51,58,61,63,68,69,70,72,77,79,80,87,89] #no 83,86
    M_red  = [11,23,24,41,52,66]
    halo   = M_blue + M_red
    therange = 'COShalos'
    if sys.argv[1] == 'halo':
        ctype = 'Mhalo'
        logmax = 12.8 #halo
        logmin = 11.5 #halo
        cblabel = r'M$_{halo}$'
    else:

        ctype = 'Mstar'
        logmax = 11.4 #star
        logmin = 10.5 #star
        cblabel = r'M$_{star}$'
elif sys.argv[1] == 'MW':
    #### MW mass Range from Michael
    halo = np.arange(25,56,1)
    therange = 'MW'
    if sys.argv[1] == 'halo':

        ctype = 'Mhalo'
        logmax = 12.2 #halo
        logmin = 11.9 #halo
        cblabel = r'M$_{halo}$'
    else:
        ctype = 'Mstar'
        logmax = 11.1 #star
        logmin = 10.2 # star
        cblabel = r'M$_{star}$'

elif sys.argv[1] == 'GMs':
    #### GM stellar mass range (10^10 - 5^10)
#    halo = [72,54,87,66,77,79,68,37,69,61,80,70,86,89,72,54,87,75,78,96,71,94,91,73,65,67,88,85,81,39,90,60,98,93,62,84,92,95,56,97]
    halo = np.arange(60,223)
    therange = 'GMstellarrange'
    if sys.argv[1] == 'halo':
        ctype = 'Mvir'
        logmax = 12.2 #halo
        logmin = 11.5 #halo
        cblabel = 'M$_{halo}$'

    else:
        ctype = 'Mstar'
        logmax = 10.7 #star
        logmin = 10 #star
        cblabel = 'M$_{star}$'
else:
    print('Invalid Option. Accepted ranges are:')
    print('COShalos, MW, or GMs')

labels = map(str, halo)
labels = list(labels)

ROMcmap = cm.get_cmap('rainbow')

# MUST HAVE TANGOS DATABASE POINTING AT ROM25 FOR THIS APRT
import tangos
#tangos.all_simulations()
#tangos.get_simulation('cosmo25').timesteps
tangos.get_halo('cosmo25/cosmo25p.768sg1bwK1BHe75.006912/halo_9')
ROM_BH_masses = []
ROM_BH_ids = []
for i in range(len(halo)):
    ROM_h = tangos.get_halo('cosmo25/cosmo25p.768sg1bwK1BHe75.006912/halo_'+labels[i])
    try:
        ROM_BHmass = ROM_h['BH_central'][0]['BH_mass']
        ROM_BH_masses.append(ROM_BHmass)
        ROM_BH_ids.append(halo[i])
        print('Halo ',halo[i],'has a Mvir =',ROM_h['Mvir'],' Mstar =',ROM_h['Mstar'],', Rvir =',ROM_h['max_radius'])
        if np.log10(ROM_h['Mstar']) < 10:
            print('M_star too LOW')
            continue

        np.savetxt('ROM_data_forhalos/ROM_BH_masses.txt',ROM_BH_masses)
        np.savetxt('ROM_data_forhalos/ROM_BH_ids.txt',ROM_BH_ids)

        ROM_Novi = np.loadtxt('/home1/nnsanche/ROMULUS_research/OVI_columndensities/COS_halos_range_data/Novi_profile_np/'+str(halo[i])+'_Novi_profile_Novi_6912.np',dtype=float)
    except:
        print('NO NOVI')
        continue
    ROM_Rbins = np.loadtxt('/home1/nnsanche/ROMULUS_research/OVI_columndensities/COS_halos_range_data/Novi_profile_np/'+str(halo[i])+'_Novi_profile_rbins_6912.np',dtype=float)


    colors = ROMcmap((np.log10(ROM_h[ctype])-logmin)/(logmax-logmin))

    SFR_at25Myr = ROM_h['SFR_encl_25Myr']/(0.6*ROM_h['Mstar'])
#    if halo[i] in M_red:#SFR_at25Myr[-1:] < 1.6*10**-11 :
#        SFR = ROM_h['inner_SFR_histogram']
#        SFR_property_object = ROM_h.get_objects('SFR_histogram')[0]
#        SFR_time_bins = SFR_property_object.x_values()
#        p.plot(SFR_time_bins, SFR)    
#        print('THIS GALAXY IS RED')        
#        plt.plot(np.log10(ROM_BHmass),np.log10(ROM_Novi[ROM_Rbins == 76.3195]),color=colors,alpha=0.5,marker='s',markersize=8)
#    else:
    plt.plot(np.log10(ROM_BHmass),np.log10(ROM_Novi[ROM_Rbins == 76.3195]),alpha=0.5,marker='o',markersize=8,color=colors)

sm = plt.cm.ScalarMappable(cmap=ROMcmap, norm=plt.Normalize(vmin=logmax, vmax=logmin))
sm._A = []

labels = ['P0','GM1','GM7','GM4']
NEW_lab = ['P0','GM1','GM2','GM3']
colors = ['DodgerBlue','SteelBlue','FireBrick','Salmon']#''Red','Salmon','Orange']
labels_noBHs = ['P0noBH','GM1noBH','GM7noBH','GM4noBH']
NEW_lab_noBHs = ['P0_noBH','GM1_noBH','GM7_noBH','GM4_noBH']
colors_noBHs = ['DodgerBlue','SteelBlue','FireBrick','Salmon']#,'Orange']
spot = ['o','o','s','s']

time = '3456'

plt.plot([0],[1],marker='o',linestyle=None,label='SF',color='Black')
plt.plot([0],[1],marker='s',linestyle=None,label='Passive',color='Black')
for k in range(len(labels)):
    print('Read in: ',labels[k])
    Novi = np.loadtxt('../profiles_plot7/'+labels[k]+'/'+labels[k]+'_Novi_'+time+'_'+catalogue+'.np')
    R = np.loadtxt('../profiles_plot7/'+labels[k]+'/'+labels[k]+'_Rbins_'+time+'_'+catalogue+'.np')
    # Read in numpy array 
    BH_times = np.loadtxt('/home1/nnsanche/PIONEER_research2/BH_studies/'+NEW_lab[k]+'_BHtime.txt')
    BH_mass  = np.loadtxt('/home1/nnsanche/PIONEER_research2/BH_studies/'+NEW_lab[k]+'_BHmass.txt')

    plt.plot(np.log10(BH_mass[(BH_times > 11.645) & (BH_times < 11.65)][0]),Novi[R == 76.],label=NEW_lab[k],color=colors[k],marker=spot[k],markersize=8,linestyle=None,markeredgecolor='Black')

plt.text(6.35,15,'z = 0.17',size=15)
plt.ylabel(r'N$_{Ovi}$ at 75 kpc',size=15)
plt.xlabel(r'M$_{BH}$/M$_{\odot}$',size=15)
plt.ylim(10.5,15.5)
plt.xlim(6,9.5)
plt.legend(fontsize=11,ncol=1,loc=4)
plt.colorbar(sm,label=r'log '+cblabel+'/M$_{\odot}$')
plt.savefig('BHmass_Noviat75kpc_'+therange+'_color'+ctype+'.pdf')
plt.show()



