# This file reads in the particles found in traced_H1gas.py, i.e. the gas particles 
# and parent gas particles of the stars in the GM1 main halo at z = 0, and determines
# where those gas particles are compared to the GM4 whole simulation

# follow this script with: makeGM1h1gasinGM4plot_temp.py: 
#           Makes pynbody plots of where the GM1 particles are in the GM4 sim & v.v.


# N. Nicole Sanchez -- June 23 2017
# Univ. of Wash.    -- Nbody Shop  
import matplotlib.pyplot as plt
import numpy as np
import pynbody

GM1h1z0iords = np.loadtxt('GM1_h1_tracedgasiord.txt')
print(GM1h1z0iords)

timesteps = np.loadtxt('timesteps.txt',dtype='str')
i = 21 # final step in sim

GM4_sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00'+timesteps[i])
# Centering on main halo of GM4
h4 = GM4_sim.halos()
pynbody.analysis.halo.center(h4[1],mode='ssc')
h4 = h4[1]
pynbody.analysis.angmom.faceon(h4)
GM4_sim.physical_units()

GM1_sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00'+timesteps[i])
# Centering on main halo of GM1
h1 = GM1_sim.halos()
pynbody.analysis.halo.center(h1[1],mode='ssc')
h1 = h1[1]
pynbody.analysis.angmom.faceon(h1)
GM1_sim.physical_units()

# np.in1d is one of the most useful things I have ever learned
# matched items in different arrays
GM1h1z0iords = np.random.choice(GM1h1z0iords,size=5000)

GM4_GM1h1z0_gasmask  = np.in1d(GM4_sim.g['iord'],GM1h1z0iords,assume_unique=False)
GM4_GM1h1z0_starmask = np.in1d(GM4_sim.s['igasorder'],GM1h1z0iords,assume_unique=False) 
print(len(GM4_sim.g['iord'][GM4_GM1h1z0_gasmask]))

GM1h1z0_gasmask  = np.in1d(GM1_sim.g['iord'],GM1h1z0iords,assume_unique=False)
GM1h1z0_starmask = np.in1d(GM1_sim.s['igasorder'],GM1h1z0iords,assume_unique=False)
print(len(GM1_sim.g['iord'][GM1h1z0_gasmask]))

#test = np.random.choice(GM4_sim.g['pos'],size=10000)
#test2 = np.random.choice(GM1_sim.g['pos'],size=10000)

GM4_gaspos = np.array(GM4_sim.g['pos'][GM4_GM1h1z0_gasmask])
GM4_gastemp = np.array(GM4_sim.g['temp'][GM4_GM1h1z0_gasmask])
GM4_gasmetal = np.array(GM4_sim.g['Metalsdot'][GM4_GM1h1z0_gasmask])
GM4_gasR = (GM4_gaspos[:,0]**2 + GM4_gaspos[:,1]**2 + GM4_gaspos[:,2]**2)**0.5

GM1_gaspos = np.array(GM1_sim.g['pos'][GM1h1z0_gasmask])
GM1_gastemp = np.array(GM1_sim.g['temp'][GM1h1z0_gasmask])
GM1_gasmetal = np.array(GM1_sim.g['Metalsdot'][GM1h1z0_gasmask])
print(GM1_sim.g['Metalsdot'])
GM1_gasR = (GM1_gaspos[:,0]**2 + GM1_gaspos[:,1]**2 + GM1_gaspos[:,2]**2)**0.5

matchgas_mask = np.in1d(GM4_sim.g['iord'][GM4_GM1h1z0_gasmask],GM1_sim.g['iord'][GM1h1z0_gasmask])
new_GM4_gastemp = np.array(GM4_gastemp[matchgas_mask])
new_GM4_gasmetal = np.array(GM4_gasmetal[matchgas_mask])
new_GM4_gaspos = np.array(GM4_gaspos[matchgas_mask])
new_GM4_gasR = (new_GM4_gaspos[:,0]**2 + new_GM4_gaspos[:,1]**2 + new_GM4_gaspos[:,2]**2)**0.5

matchgas_mask1 = np.in1d(GM1_sim.g['iord'][GM1h1z0_gasmask],GM4_sim.g['iord'][GM4_GM1h1z0_gasmask])
new_GM1_gastemp = np.array(GM1_gastemp[matchgas_mask1])
new_GM1_gasmetal = np.array(GM1_gasmetal[matchgas_mask1])
new_GM1_gaspos = np.array(GM1_gaspos[matchgas_mask1])
new_GM1_gasR = (new_GM1_gaspos[:,0]**2 + new_GM1_gaspos[:,1]**2 + new_GM1_gaspos[:,2]**2)**0.5

#plt.plot(new_GM4_gasR,new_GM1_gasR,marker='o',linestyle='None',markersize=1)
#plt.xlabel('GM4 Gas Distance from Center [kpc]')
#plt.ylabel('GM1 Gas Distance from Center [kpc]')
#plt.savefig('compareR_GM1_GM4.pdf')
#plt.show()

#plt.plot(new_GM4_gaspos[:,0],new_GM4_gaspos[:,1],linestyle='None',markersize=1,marker='o')
np.savetxt('GM4_GM1matchedgas_posxyz.txt',new_GM4_gaspos)
np.savetxt('GM4_GM1matchedgas_temp.txt',new_GM4_gastemp)
np.savetxt('GM4_GM1matchedgas_metal.txt',new_GM4_gasmetal)
np.savetxt('GM1matchedgas_posxyz.txt',new_GM1_gaspos)
np.savetxt('GM1matchedgas_temp.txt',new_GM1_gastemp)
np.savetxt('GM1matchedgas_metal.txt',new_GM1_gasmetal)



