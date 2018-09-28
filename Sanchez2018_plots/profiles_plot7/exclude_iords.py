import pynbody
import numpy as np
import sys

time = str(sys.argv[1])
halo = int(sys.argv[2])

sim = pynbody.load('/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.00'+time)
#sim = pynbody.load('/nobackup/nnsanche/NO_BHs/pioneer50h243GM7.1536gst1bwK1/pioneer50h243GM7.1536gst1bwK1.00'+time)
#sim = pynbody.load('/nobackup/nnsanche/NO_BHs/pioneer50h243GM4.1536gst1bwK1/pioneer50h243GM4.1536gst1bwK1.00'+time)
h   = sim.halos()
sat = h[halo]
sat_iords = sat.g['iord']
print('Number of particles in the satellite',str(halo),' at timestep',time,':',len(sat_iords))
print('Saving in text file')
np.savetxt('GM7_h'+str(halo)+'iords_'+time+'.txt',sat_iords)
