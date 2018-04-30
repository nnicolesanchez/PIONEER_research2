# This script picks which two particles I want to track across all timsteps

# See script: track2particles_GM1_GM4.py for:
#      Both particles:
#           x, y, z, v, rho, temp, metals
#           ** When stars set rho & temp = 0 
#           ** Multiply xyz * v by expansion factor

# These are two particles which for stars in GM1 but
# remain only gas (form no stars) in GM4


# N. Nicole Sanchez -- July 1 2017
# Univ. of Wash.    -- Nbody Shop  
import matplotlib.pyplot as plt
import numpy as np
import pynbody

#GM1h1z0iords = np.loadtxt('GM1_h1_tracedgasiord.txt')
GM1h1z0igasorder = np.loadtxt('GM1_h1_tracedigasorder.txt') 

timesteps = np.loadtxt('timesteps.txt',dtype='str')
i = 21 # final step in sim

GM4_sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00'+timesteps[i])
h4 = GM4_sim.halos()
h4 = h4[1]
pynbody.analysis.halo.center(h4,mode='ssc')
pynbody.analysis.angmom.faceon(h4)
GM4_sim.physical_units()

GM1_sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00'+timesteps[i])
h1 = GM1_sim.halos()
h1 = h1[1]
pynbody.analysis.halo.center(h1,mode='ssc')
pynbody.analysis.angmom.faceon(h1)
GM1_sim.physical_units()

# Get two star particles from GM1 that form in the last 2-3 Gyr (once GM4 quenches)
# i.e. formed after ~7 Gyr
GM1_sim.s['tform'].in_units('Gyr')

stars1 = GM1_sim.s
youngstars1 = h1.s[h1.s['tform'].in_units('Gyr') > 7]
youngstarsigasords1 =  youngstars1['igasorder']

gas4 = GM4_sim.g
gas4_parents_mask = np.in1d(GM4_sim.g['iord'],youngstarsigasords1)
gas4_parentsof1 = GM4_sim[gas4_parents_mask]
most_mass = gas4_parentsof1[np.where(gas4_parentsof1.g['mass'].in_units('Msol') > 5*10**5)]

print('Most massive gas particles in GM4 that match parent gas particles of stars in GM1')
print('I.e. most likely to have NOT formed any stars in GM4 by z=0')
print(most_mass.g['iord'])

particles = most_mass.g['iord']
#particles = particles[[1,3,4,6,7,8,9,11]]
#print(particles)

for i in range(len(particles)):
    print('Particle ',particles[i],'forms this many stars in GM4:')
    print('(The following array should be empty for me to track it)')
    print(GM4_sim.s['iord'][GM4_sim.s['igasorder'] == particles[i]])
    

# Gonna trace these particles through time!
# Only saving the ones that had empty star particle arrays in GM4
np.savetxt('tracetheseparticles.txt',particles[[1,3,4,6,7,8,9,11]])


# Choose two random particles in following script:  track2particles_GM1_GM4.py
# Like the following examples

# Particle 1: 881021  
# Definitely doesn't make stars in GM4
# Turns into 2 stars in GM1
#GM1_sim.s['iord'][GM1_sim.s['igasorder'] == particles[1]]

# Particle 2: 594778 
# Definitely doesn't make stars in GM4       
# Turns into 4 stars in GM1
#GM1_sim.s['iord'][GM1_sim.s['igasorder'] == particles[0]]
