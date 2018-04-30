# This script creates a phase diagram (temp vs density)
# for the CGM of a ChaNGa Nbody Simulation

#     - Outputs:
#          plots of phase diagram


# N. Nicole Sanchez -- July 2 2017
# Univ. of Wash.    -- Nbody Shop
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np
import pynbody
import sys

sim = ['/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.00']
labels = ['P0','GM1','GM4','GM5','GM6','GM7']



if len(sys.argv) == 1:
    print('No galaxy selected. Current options: P0, GM1, GM4, GM5, GM6')
    print('Syntax: "CGMphasediagram.py GM1"')
    quit()
else:
    if (str(sys.argv[1]) == 'P0'):
        k = 0 
    elif (str(sys.argv[1]) == 'GM1'):
        k = 1
    elif (str(sys.argv[1]) == 'GM4'):
        k = 2
    elif (str(sys.argv[1]) == 'GM5'):
        k = 3
    elif (str(sys.argv[1]) == 'GM6'):
        k = 4
    elif (str(sys.argv[1]) == 'GM7'):
        k = 5
    else :
        print('Not a valid option. Current options: P0, GM1, GM4, GM5, GM6, GM7')
        print('Syntax: "CGMphasediagram.py GM1"')
        quit()    

    ts = np.loadtxt('../'+labels[k]+'/timesteps.txt',dtype=str)
    t = len(ts)-1
    if k == 5:
        t = len(ts)-2
    sim = pynbody.load(sim[k]+ts[t])
    name = str(sys.argv[1])
    print(name+' simulation at z = ','%.2f' % sim.properties['z'] )

sim.properties
sim.properties['time'].in_units('Gyr')
sim.loadable_keys()
sim.physical_units()
h = sim.halos()
h1 = h[1]
pynbody.analysis.angmom.faceon(h1)

# Constants
m_H = 1.6733 * 10**-24 #g

# Want to isolate CGM  
# Isolate and remove disk stars within radius 0-10 kpc & vertically 4 kpc 
r_max = 10  # kpc
z_max = 10   # kpc

Rg_d = ((h1.g['x'].in_units('kpc'))**2. + (h1.g['y'].in_units('kpc'))**2.)**(0.5)
disk_gas_xymax =  (Rg_d < r_max)
disk_gas_zmax  = (h1.g['z'].in_units('kpc') < z_max) & (h1.g['z'].in_units('kpc') > -z_max)

disk_gas_mask = disk_gas_xymax & disk_gas_zmax
disk_gas = h1.g[disk_gas_xymax & disk_gas_zmax]
CGM_gas  = h1.g[~disk_gas_mask]

x = np.log10(CGM_gas['rho'].in_units('g cm**-3')/m_H)
y = np.log10(CGM_gas['temp'])
z = np.log10(CGM_gas['mass'])

# Specifically isolate the cool/dense tail
x = x[y < 4.4]
z = z[y < 4.4]
y = y[y < 4.4]

pynbody.plot.sph.faceon_image(CGM_gas, qty='rho', width='200 kpc', resolution=500,units='g cm^-3')
plt.show()

pynbody.plot.sph.sideon_image(CGM_gas, qty='rho', width='200 kpc', resolution=500,units='g cm^-3')
plt.show()

pynbody.plot.sph.sideon_image(CGM_gas[CGM_gas['temp'] < 10**4.4], qty='rho', width='200 kpc', resolution=500,units='g cm^-3')
plt.show()

pynbody.plot.sph.faceon_image(CGM_gas[CGM_gas['temp'] < 10**4.4], qty='rho', width='200 kpc', resolution=500,units='g cm^-3')
plt.show()

#fig = plt.figure(figsize=(7, 5))
#plt.hist2d(x,y,(100,100),weights=z,cmap=cm.jet,norm=mpl.colors.LogNorm())
#plt.ylabel(r'Log$_{10}$ T ('+str(CGM_gas['temp'].units)+')')
#plt.xlabel(r'Log$_{10}$ n$_H$ (cm$^{-3}$)')
#plt.colorbar(label=(r'Log M (M$_{\odot}$)'))
#plt.text(-5.5,6.7,name, color='midnightblue',size=12)
#plt.text(0.7,6.7,'z = 0',color='midnightblue',size=12)
#plt.xlim(-6,2)
#plt.ylim(3.5,7)
#plt.savefig(name+'_phasediagram.pdf')
#plt.show()

