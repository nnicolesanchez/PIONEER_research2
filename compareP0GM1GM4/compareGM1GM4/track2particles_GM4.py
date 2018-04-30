# This script reads in tracetheseparticles.txt which contains iords 
# for gas particles that are 
#       - stars in GM1 by z=0
#       - still gas and never form stars in GM4
#       - form stars in GM1 after 7 Gyr (i.e. after quenching occured in GM4)

# Creates trace+########_GM#.csv file with the following info for 
# 1 'tracetheseparticles.txt' particles across
# all timesteps!
#       - mass, pos(xyz), vel(xyz), rho, temp, metal fraction


# Need to run: pick2particles_GM1star_GM4gas.py first
# To create: tracetheseparticles.txt


# N. Nicole Sanchez -- July 1 2017
# Univ. of Wash.    -- Nbody Shop
import matplotlib.pyplot as plt
import numpy as np
import pynbody
import csv

# Change this to look at other particles 
# len(particles) = 7
N = 0

particles = np.loadtxt('tracetheseparticles.txt')
print(len(particles))
p1 = particles[N]
print('This is the gas iord of the particle I am tracking:', p1)

timesteps = np.loadtxt('timestepsGM4_GM1_numbers.txt',dtype='str')
#i = 36 # final step in sim


with open('GM4_'+str(int(p1))+'.csv', 'w') as csvfile:
    fieldnames = ['time','redshift','iords', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'rho', 'temp', 'metals']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=',')

    writer.writeheader()
    
    for i in range(len(timesteps)):
        print(timesteps[i])
        GM4_sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.'+timesteps[i])
        h4 = GM4_sim.halos()
        if (len(h4) == 0):
            continue
        #h4 = h4[1]
        pynbody.analysis.halo.center(h4[1],mode='ssc')
        pynbody.analysis.angmom.faceon(h4[1])
        GM4_sim.physical_units()
        
        p1_deets = GM4_sim.g[GM4_sim.g['iord'] == p1]
        
        p1_iord = p1_deets['iord']
        p1_mass = p1_deets['mass']
        p1_pos  = p1_deets['pos'][0]*GM4_sim.properties['a']
        p1_vel  = p1_deets['vel'][0]*GM4_sim.properties['a']
        p1_rho  = p1_deets['rho']
        p1_temp = p1_deets['temp']
        p1_met  = p1_deets['metals']
        print(GM4_sim.properties['a'])
        print(p1_iord,p1_mass,p1_pos,p1_vel,p1_rho,p1_temp,p1_met)

        if (i == 0):
            writer.writerow({'time': str('Gyr'), 'redshift': str('NoUnit()'), 'iords': str(p1_iord.units), 'mass': str(p1_mass.units), 'x': str(p1_deets['pos'].units)+'*a', 'y': str(p1_deets['pos'].units)+'*a', 'z': str(p1_deets['pos'].units)+'*a', 'vx': str(p1_deets['vel'].units)+'*a', 'vy': str(p1_deets['vel'].units)+'*a', 'vz': str(p1_deets['vel'].units)+'*a', 'rho': str(p1_rho.units),'temp': str(p1_temp.units), 'metals': str(p1_met.units)})

            writer.writerow({'time': str("%.3f" % GM4_sim.properties['time'].in_units('Gyr')), 'redshift': str("%.3f" % ((1/GM4_sim.properties['a']) - 1)), 'iords': str(p1_iord[0]), 'mass': str("%.3f" % p1_mass[0]), 'x': str("%.3f" % p1_pos[0]), 'y': str("%.3f" % p1_pos[1]), 'z': str("%.3f" % p1_pos[2]), 'vx': str("%.3f" % p1_vel[0]), 'vy': str("%.3f" % p1_vel[1]), 'vz': str("%.3f" % p1_vel[2]), 'rho': str("%.3f" % p1_rho[0]), 'temp': str("%.3f" % p1_temp[0]), 'metals': str("%.3f" % p1_met[0])})

        else: 
            writer.writerow({'time': str("%.3f" % GM4_sim.properties['time'].in_units('Gyr')), 'redshift': str("%.3f" % ((1/GM4_sim.properties['a']) - 1)), 'iords': str(p1_iord[0]), 'mass': str("%.3f" % p1_mass[0]), 'x': str("%.3f" % p1_pos[0]), 'y': str("%.3f" % p1_pos[1]), 'z': str("%.3f" % p1_pos[2]), 'vx': str("%.3f" % p1_vel[0]), 'vy': str("%.3f" % p1_vel[1]), 'vz': str("%.3f" % p1_vel[2]), 'rho': str("%.3f" % p1_rho[0]), 'temp': str("%.3f" % p1_temp[0]), 'metals': str("%.3f" % p1_met[0])})



