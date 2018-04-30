# This script plots the star formation history for Patient 0 
import matplotlib.pyplot as plt
import numpy as np
import pynbody

timesteps = np.loadtxt('timesteps.txt',dtype='str')
plt.plot([0,0.],[0,0.01],color='SteelBlue',label='GM1')
plt.plot([0,0],[0,0.01],color='FireBrick',label='GM4')

#for j in range(1):#len(timesteps)):
i = 21
file = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00'+timesteps[i])

h = file.halos()
pynbody.analysis.halo.center(h[1],mode='ssc')
h1 = h[1]

pynbody.analysis.angmom.faceon(h1)
file.physical_units()

sfhist, bins = pynbody.plot.stars.sfh(file,filename='GM4_sfh_xlim'+timesteps[i]+'.pdf',massform=True,color='FireBrick',legend=True,trange=[0,14])
#pynbody.plot.stars.sfh(h1,filename='GM4_MH_sfh_'+timesteps[i]+'.pdf',massform=False)
    
file2 = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00'+timesteps[i])

h_2 = file2.halos()
pynbody.analysis.halo.center(h_2[1],mode='ssc')
h1_2 = h_2[1]

pynbody.analysis.angmom.faceon(h1_2)
file.physical_units()
    #pynbody.plot.stars.sfh(h1,filename='GM1_MH_sfh_'+timesteps[i]+'.pdf',massform=False)

plt.ylim(0,12)
#    plt.xlim(0,14)
sfhist2, bins2 = pynbody.plot.stars.sfh(file2,filename='GM1_4_sfh_xlim'+timesteps[i]+'.pdf',massform=True,legend=True,color='SteelBlue',trange=[0,14])
#plt.xlim(0,15)
plt.legend(loc=2)
plt.savefig('GM1_4_sfh_'+timesteps[i]+'.pdf')
plt.show()
print(timesteps[i])

np.savetxt('GM4_sfhistory_bins.txt', np.transpose([sfhist, bins[:-1]]))
np.savetxt('GM1_sfhistory_bins.txt', np.transpose([sfhist2, bins2[:-1]]))
