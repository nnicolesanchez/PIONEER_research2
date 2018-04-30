# This script plots the star formation history for Patient 0 
import matplotlib.pyplot as plt
import numpy as np
import pynbody

timesteps = np.loadtxt('timesteps.txt',dtype='str')


for j in range(1):#len(timesteps)):
    i = 21
    file = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00'+timesteps[i])

    h = file.halos()
    pynbody.analysis.halo.center(h[1],mode='ssc')
    h1 = h[1]

    pynbody.analysis.angmom.faceon(h1)
    file.physical_units()

    pynbody.plot.stars.sfh(file,filename='GM1_sfh_xlim'+timesteps[i]+'.pdf',massform=True)
    #pynbody.plot.stars.sfh(h1,filename='GM1_MH_sfh_'+timesteps[i]+'.pdf',massform=False)
    
    #plt.xlim(0,14)
    plt.show()
    print(timesteps[i])

