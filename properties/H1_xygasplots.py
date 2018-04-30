# This script plots face-on and edge-on plots of ChaNGa Nbody
# simulations in gas density and overplots the virial radius


# N. Nicole Sanchez -- July 5, 2017
# U. W. Seattle     -- Nbody Shop
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import numpy as np
import pynbody
import sys

plt.rc('font', size=12, family='serif', style='normal', variant='normal', stretch='normal', weight='normal')
plt.rc('xtick', labelsize=12)
plt.rc('xtick.major', size=6, width=1)
plt.rc('lines', lw=2)
plt.rc('axes', lw=1, labelsize=12)

if len(sys.argv) == 1:
    print('No galaxy selected. Current options: P0, GM1, GM4, GM5, GM6, GM7')
    print('Syntax: "GM_xygasplots.py GM1"')
    quit()
else:
    if (str(sys.argv[1]) == 'P0'):
        sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.004096')
    elif (str(sys.argv[1]) == 'GM1'):
        sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.004096')
    elif (str(sys.argv[1]) == 'GM4'):
        sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/OLD/pioneer50h243GM4.1536gst1bwK1BH.004096')
    elif (str(sys.argv[1]) == 'GM5'):
        sim = pynbody.load('/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.004096')
    elif (str(sys.argv[1]) == 'GM6'):
        sim = pynbody.load('/nobackupp8/fgoverna/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.004096')
    elif (str(sys.argv[1]) == 'GM7'):
        sim = pynbody.load('/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.004096')

    else :
        print('Not a valid option. Current options: P0, GM1, GM4, GM5, GM6, GM7')
        print('Syntax: "GM_xygasplots.py GM1"')
        quit()    

    name = str(sys.argv[1])
    print(name+' simulation at z = ','%.2f' % sim.properties['z'] )
    R_vir = np.loadtxt('../'+name+'/'+name+'_mainhalo.stat',dtype=float,skiprows=1,usecols=6,unpack=True)

if (str(sys.argv[1]) == 'GM7'):
    R_vir = h1.properties['Rvir'] 

print('Virial radius:',float(R_vir))
R_vir = float(R_vir)

h = sim.halos()
h1 = h[1]
pynbody.analysis.halo.center(h1,mode='ssc')
pynbody.analysis.angmom.faceon(h1)
sim.physical_units()

m_h = 1.6733 * 10**-24 # g
sim.g['rho_in_nhcm'] = sim.g['rho'].in_units('g cm**-3') / m_h


# Plotting face on
pynbody.plot.sph.image(sim.g, qty='rho_in_nhcm' , width=2*270, cmap='jet', show_cbar=False, vmin=10**-6, vmax=10**1)
fig = plt.gcf()
ax = fig.gca()


circle2 = plt.Circle((0.0, 0.0), R_vir, color='white',fill=False,linestyle='--')
ax.add_artist(circle2)
cbar = plt.colorbar()
cbar.set_label(r'n$_H$ [cm$^{-3}$]')
plt.text(-245,220,name, color='White')
plt.text(185,-240,'z = 0',color='White')
#plt.xlim(275,275)
#plt.ylim(275,275)
plt.savefig(name+'_gasdensity.pdf')
plt.show()
plt.clf()

quit()
# Plotting temperature map
pynbody.plot.sph.image(sim.g, qty='temp', width=2*R_vir, cmap='jet', show_cbar=False,vmin=10**3, vmax=10**7)
fig = plt.gcf()
ax = fig.gca()

circle2 = plt.Circle((0.0, 0.0), R_vir, color='Black',fill=False,linestyle='--')
ax.add_artist(circle2)
cbar = plt.colorbar()
cbar.set_label(r'T [K]')
plt.text(-245,220,name, color='Black')
plt.text(185,-240,'z = 0',color='Black')
plt.savefig(name+'_gastemp.pdf')
plt.show()

# Plotting metal mass fraction
pynbody.plot.sph.image(sim.g, qty='metals', width=2*R_vir, cmap='jet', show_cbar=False,vmin=10**-6, vmax=5*10**-1)
fig = plt.gcf()
ax = fig.gca()

circle2 = plt.Circle((0.0, 0.0), R_vir, color='Black',fill=False,linestyle='--')
ax.add_artist(circle2)
cbar = plt.colorbar()
cbar.set_label(r'metal mass fraction')
plt.text(-245,220,name, color='Black')
plt.text(185,-240,'z = 0',color='Black')
plt.savefig(name+'_gasmetalfrac.pdf')
plt.show()
