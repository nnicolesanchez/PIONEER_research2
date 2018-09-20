#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

#catalogue = 'grp'
catalogue = 'AHF'

### IDs for galaxies with stellar masses within COS-Halos range
M_blue = [9,14,15,16,19,20,21,22,27,28,29,30,31,32,34,38,41,42,43,44,46,47,48,51,55,56,60,61,62,64,65,66,71,73,75,77,79,81,82,92] #lost 12,74 at 6912 
M_red  = [10,11,17,18,23,24,26,33,35,36,37,45,49,50,58,70] 
halo = M_blue + M_red
labels = map(str, halo)
labels = list(labels)

# MUST HAVE TANGOS DATABASE POINTING AT ROM25 FOR THIS APRT
import tangos
#tangos.all_simulations()
#tangos.get_simulation('cosmo25').timesteps
tangos.get_halo('cosmo25/cosmo25p.768sg1bwK1BHe75.006912/halo_9')
ROM_BH_masses = []
ROM_BH_ids = []
for i in range(len(halo)):
    ROM_h = tangos.get_halo('cosmo25/cosmo25p.768sg1bwK1BHe75.006912/halo_'+labels[i])
    ROM_BHmass = ROM_h['BH_central'][0]['BH_mass']
    ROM_BH_masses.append(ROM_BHmass)
    ROM_BH_ids.append(halo[i])
    print('Halo ',halo[i],'has a Mvir =',ROM_h['Mvir'],' Mstar =',ROM_h['Mstar'],', Rvir =',ROM_h['max_radius'])

    np.savetxt('ROM_data_forhalos/ROM_BH_masses.txt',ROM_BH_masses)
    np.savetxt('ROM_data_forhalos/ROM_BH_ids.txt',ROM_BH_ids)

    ROM_Novi = np.loadtxt('/home1/nnsanche/ROMULUS_research/OVI_columndensities/COS_halos_range_data/Novi_profile_np/'+str(halo[i])+'_Novi_profile_Novi_6912.np',dtype=float)
    ROM_Rbins = np.loadtxt('/home1/nnsanche/ROMULUS_research/OVI_columndensities/COS_halos_range_data/Novi_profile_np/'+str(halo[i])+'_Novi_profile_rbins_6912.np',dtype=float)

#    print(ROM_Rbins)
    SFR_at25Myr = ROM_h['SFR_encl_25Myr']/(0.6*ROM_h['Mstar'])
    if halo[i] in M_red:#SFR_at25Myr[-1:] < 1.6*10**-11 :
#        SFR = ROM_h['inner_SFR_histogram']
#        SFR_property_object = ROM_h.get_objects('SFR_histogram')[0]
#        SFR_time_bins = SFR_property_object.x_values()
#        p.plot(SFR_time_bins, SFR)    
        print('THIS GALAXY IS RED')        
        plt.plot(np.log10(ROM_BHmass),np.log10(ROM_Novi[ROM_Rbins == 76.3195]),color='DarkGrey',alpha=0.5,marker='s',markersize=8)
    else:
        plt.plot(np.log10(ROM_BHmass),np.log10(ROM_Novi[ROM_Rbins == 76.3195]),color='DarkGrey',alpha=0.5,marker='o',markersize=8)


labels = ['P0','GM1','GM7','GM4']
NEW_lab = ['P0','GM1','GM2','GM3']
colors = ['DodgerBlue','SteelBlue','FireBrick','Salmon']#''Red','Salmon','Orange']
labels_noBHs = ['P0noBH','GM1noBH','GM7noBH','GM4noBH']
NEW_lab_noBHs = ['P0_noBH','GM1_noBH','GM7_noBH','GM4_noBH']
colors_noBHs = ['DodgerBlue','SteelBlue','FireBrick','Salmon']#,'Orange']
spot = ['o','o','s','s']

time = '3456'

plt.plot([0],[1],marker='o',linestyle=None,label='Star Forming',color='Black')
plt.plot([0],[1],marker='s',linestyle=None,label='Passive',color='Black')
for k in range(len(labels)):
    print('Read in: ',labels[k])
    Novi = np.loadtxt('../profiles_plot7/'+labels[k]+'/'+labels[k]+'_Novi_'+time+'_'+catalogue+'.np')
    R = np.loadtxt('../profiles_plot7/'+labels[k]+'/'+labels[k]+'_Rbins_'+time+'_'+catalogue+'.np')
    # Read in numpy array 
    BH_times = np.loadtxt('/home1/nnsanche/PIONEER_research2/BH_studies/'+NEW_lab[k]+'_BHtime.txt')
    BH_mass  = np.loadtxt('/home1/nnsanche/PIONEER_research2/BH_studies/'+NEW_lab[k]+'_BHmass.txt')

    plt.plot(np.log10(BH_mass[(BH_times > 11.645) & (BH_times < 11.65)][0]),Novi[R == 76.],label=NEW_lab[k],color=colors[k],marker=spot[k],markersize=8,linestyle=None)

plt.text(6.75,15.25,'z = 0.17',size=15)
plt.ylabel(r'N$_{Ovi}$ at 75 kpc',size=15)
plt.xlabel(r'M$_{BH}$/M$_{\odot}$',size=15)
plt.ylim(13.,15.5)
plt.xlim(6.6,9.5)
plt.legend(fontsize=12)#,loc=4)
plt.savefig('BHmass_Noviat75kpc.pdf')
plt.show()



