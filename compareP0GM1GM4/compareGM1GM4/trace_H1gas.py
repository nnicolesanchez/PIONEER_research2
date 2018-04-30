# This script determines which gas and star particles end up in the main
# halo at z = 0 of GM1 (the active one with slightly smaller satellite 
# than Patient 0) and the parent gas particles of all the stars too.

# Prints out all iords for gas & parent gas particles in GM1

import matplotlib.pyplot as plt
import numpy as np
import pynbody

timesteps = np.loadtxt('timesteps.txt',dtype='str')

i = 21
sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00'+timesteps[i])
#sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00'+timesteps[i])

sim.properties
sim.loadable_keys()

sim.physical_units

h = sim.halos()
h1 = h[1]

gasatz_0 = h1.g['iord']
starsatz_0_gasparents = h1.s['igasorder']

set(h1.g['iord']) & set(h1.s['igasorder'])
overlap = set(h1.g['iord']) & set(h1.s['igasorder'])

print('Number of gas particles in H1 at z = 0: ',len(h1.g['iord']))
print('Number of star particles in H1 at z = 0: ',len(h1.s['iord']))
print('Number of gas parents for stars in H1 at z = 0: ',len(h1.s['igasorder']))

h1gas = h1.g['iord']
h1stargas = h1.s['igasorder']

h1tracedgas = np.unique(h1gas)
h1tracedgas = [int(i) for i in h1tracedgas]  
np.savetxt('GM1_h1_tracedgasiord.txt',h1tracedgas,fmt='%07d')

h1tracedstar = np.unique(h1stargas)
h1tracedstar = [int(i) for i in h1tracedstar]   
np.savetxt('GM1_h1_tracedigasorder.txt',h1tracedstar,fmt='%07d')

h1allgas_history = np.concatenate((h1gas,h1stargas),axis=0)
print('Total number of gas particles referenced by gas iord and star parents: ',len(h1allgas_history))
print('Total number of unique gas particles + star parents in H1: ',len(np.unique(h1allgas_history)))

#h1tracedgas = np.unique(h1allgas_history)
#h1tracedgas = [int(i) for i in h1tracedgas]
#np.savetxt('GM1_h1_tracedgasiord.txt',h1tracedgas,fmt='%07d')
