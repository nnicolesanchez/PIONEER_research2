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

plt.rc('font', size=12, family='serif', style='normal', variant='normal', stretch='normal', weight='normal')
plt.rc('xtick', labelsize=12)
plt.rc('xtick.major', size=6, width=1)
plt.rc('lines', lw=2)
plt.rc('axes', lw=1, labelsize=12)

if len(sys.argv) == 1:
    print('No galaxy selected. Current options: P0, GM1, GM4, GM5, GM6')
    print('Syntax: "GM_xygasplots.py GM1"')
    quit()
else:
    if (str(sys.argv[1]) == 'P0'):
        sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.004096')
    elif (str(sys.argv[1]) == 'GM1'):
        sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.004096')
    elif (str(sys.argv[1]) == 'GM4'):
        sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.004096')    
    elif (str(sys.argv[1]) == 'GM5'):
        sim = pynbody.load('/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.004096')
    elif (str(sys.argv[1]) == 'GM6'):
        sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.004096')
    elif (str(sys.argv[1]) == 'GM7'):
        sim = pynbody.load('/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.004096') 

    else :
        print('Not a valid option. Current options: P0, GM1, GM4')
        print('Syntax: "GM_xygasplots.py GM1"')
        quit()    

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
z_max = 4   # kpc

Rg_d = ((h1.g['x'].in_units('kpc'))**2. + (h1.g['y'].in_units('kpc'))**2.)**(0.5)
disk_gas_xymax =  (Rg_d < r_max)
disk_gas_zmax  = (h1.g['z'].in_units('kpc') < z_max) & (h1.g['z'].in_units('kpc') > -z_max)

disk_gas_mask = disk_gas_xymax & disk_gas_zmax
disk_gas = h1.g[disk_gas_xymax & disk_gas_zmax]
CGM_gas  = h1.g[~disk_gas_mask]

x = np.log10(CGM_gas['rho'].in_units('g cm**-3')/m_H)
y = np.log10(CGM_gas['temp'])
z = np.log10(CGM_gas['mass'])
#print(CGM_gas['mass'].units)

fig = plt.figure(figsize=(7, 5))
plt.hist2d(x,y,(100,100),weights=z,cmap=cm.jet,norm=mpl.colors.LogNorm())
plt.ylabel(r'Log$_{10}$ T ('+str(CGM_gas['temp'].units)+')')
plt.xlabel(r'Log$_{10}$ n$_H$ (cm$^{-3}$)')
plt.colorbar(label=(r'Log M (M$_{\odot}$)'))
plt.text(-5.5,6.7,name, color='midnightblue',size=12)
plt.text(0.7,6.7,'z = 0',color='midnightblue',size=12)
plt.xlim(-6,2)
plt.ylim(3.5,7)
plt.savefig(name+'_phasediagram.pdf')
plt.show()
