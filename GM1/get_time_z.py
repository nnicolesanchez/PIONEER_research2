import pynbody 
import numpy as np

ts = np.loadtxt('timesteps.txt',dtype='str')

z = []
gyr = []
for i in range(len(ts)):
    sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00'+ts[i])

    h = sim.halos()
    h1 = h[1]
    pynbody.analysis.halo.center(h1,mode='ssc')
    pynbody.analysis.angmom.faceon(h1)
    sim.physical_units()

    print('Timestep:',ts[i],'Redshift:',((1/sim.properties['a']) - 1),'Time:',sim.properties['time'].in_units('Gyr'))
    z.append(((1/sim.properties['a']) - 1))
    gyr.append(sim.properties['time'].in_units('Gyr'))

    np.savetxt('times_gyr.txt',gyr)
    np.savetxt('redshifts.txt',z)
