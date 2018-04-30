# Calculates the virial temperature for the main halo across time
# using the main halo information from amiga.stat
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


sims = ['/nobackupp8/fgoverna/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00','/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.00','/nobackupp8/fgoverna/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.00']
labels = ['P0','GM4']    #'GM1', 'GM4', 'GM5', 'GM6']
#'/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00',

# only 0 & 3 

s = 0



omega_M = 0.3086
omega_lam = 0.6914
for i in range(len(labels)):
    P0_Tvir = []
    P0_ts  = np.loadtxt('/u/nnsanche/PIONEER_research/'+str(labels[i])+'/timesteps.txt',dtype='str')
    for j in range(len(P0_ts)):
        P0_1 = pd.read_csv(str(sims[i])+str(P0_ts[j])+'.amiga.stat',delim_whitespace=True,nrows=3,usecols=[0,1,2,3,4,5,6])


        P0_z   = np.loadtxt('/u/nnsanche/PIONEER_research/'+str(labels[i])+'/redshifts.txt') 
        P0_Gyr = np.loadtxt('/u/nnsanche/PIONEER_research/'+str(labels[i])+'/time_Gyr.txt')


        print('Timestep:',P0_ts[j])
        M_200 = P0_1['Mvir(M_sol)'][0]
        print('M_vir:',M_200)
        T_vir = 10**5.69* (M_200/(10**12))**(2./3.) * (omega_M*(1 + P0_z[j])**3. + omega_lam)**(2./3.)
        print('T_vir: ',T_vir,'K at z = ',P0_z[j])
        P0_Tvir.append(T_vir)

    plt.plot(P0_Gyr,np.log10(P0_Tvir),label=labels[s])

plt.ylabel(r'T$_{vir}$ [log(K)]')
plt.xlabel('Time [Gyr]')
plt.legend()
plt.show()

