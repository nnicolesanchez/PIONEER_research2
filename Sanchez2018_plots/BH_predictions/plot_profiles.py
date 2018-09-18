#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

#catalogue = 'grp'
catalogue = 'AHF'

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

plt.text(7.61,14.71,'z = 0.17',size=15)
plt.ylabel(r'N$_{Ovi}$ at 75 kpc',size=15)
plt.xlabel(r'M$_{BH}$/M$_{\odot}$',size=15)
plt.ylim(14.4,14.75)
plt.xlim(7.6,7.85)
plt.legend(fontsize=15,loc=4)
plt.savefig('BHmass_Noviat75kpc.pdf')
plt.show()



