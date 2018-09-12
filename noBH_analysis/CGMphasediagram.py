# This script creates a phase diagram (temp vs density)
# for the CGM of a ChaNGa Nbody Simulation

#     - Inputs:
#          GM galaxy you want to analyze
#          timestep at which you want to analyze

#     - Outputs:
#          plots of phase diagram with colorbars in
#                   mass, metallicity, rkpc


#1;95;0c N. Nicole Sanchez -- July 2 2017
# Univ. of Wash.    -- Nbody Shop
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np
import pynbody
import sys


if len(sys.argv) == 3:
    ts = str(sys.argv[2])
else:
    print('No timestep specified. Running at z = 0.17 because reasons.')
    ts = '3456' # z= 0.17

if len(sys.argv) == 1:
    print('No galaxy selected. Current options: P0, GM1, GM7 (aka GM2), GM4 (aka GM3)')
    print('Includes galaxies without BH physics: Call "P0_noBH"')
    print('To specify timestep include 4 digit snapshot number as second input')
    print('Syntax: "CGMphasediagram.py GM1_noBH" 4096')
    quit()
else:
    if (str(sys.argv[1]) == 'P0'):
        sim = pynbody.load('/nobackupp2/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00'+ts)
    elif (str(sys.argv[1]) == 'GM1'):
        sim = pynbody.load('/nobackupp2/nnsanche/pioneer50h243GM1.1536gst1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00'+ts)
    elif (str(sys.argv[1]) == 'GM7'):
        sim = pynbody.load('/nobackupp2/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.00'+ts) 
    elif (str(sys.argv[1]) == 'GM4'):
        sim = pynbody.load('/nobackupp2/nnsanche/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00'+ts)    
    elif (str(sys.argv[1]) == 'P0noBH'):
        sim = pynbody.load('/nobackupp2/nnsanche/NO_BHs/pioneer50h243.1536gst1bwK1/pioneer50h243.1536gst1bwK1.00'+ts)
    elif (str(sys.argv[1]) == 'GM1noBH'):
        sim = pynbody.load('/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM1.1536gst1bwK1/pioneer50h243GM1.1536gst1bwK1.00'+ts)
    elif (str(sys.argv[1]) == 'GM7noBH'):
        sim = pynbody.load('/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM7.1536gst1bwK1/pioneer50h243GM7.1536gst1bwK1.00'+ts) 
    elif (str(sys.argv[1]) == 'GM4noBH'):
        sim = pynbody.load('/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM4.1536gst1bwK1/pioneer50h243GM4.1536gst1bwK1.00'+ts)    

    else :
        print('Not a valid option. Current options: P0, GM1, GM7 (aka GM2), GM4 (aka GM3)')
        print('Includes galaxies without BH physics: Call "P0_noBH"')
        print('Syntax: "CGMphasediagram.py GM1"')
        quit()    

    if str(sys.argv[1]) == 'GM7':
        name = 'GM2'
    elif str(sys.argv[1]) == 'GM4':
        name = 'GM3'
    elif str(sys.argv[1]) == 'GM7noBH':
        name = 'GM2noBH'
    elif str(sys.argv[1]) == 'GM4noBH':
        name = 'GM3noBH'
    else:
        name = str(sys.argv[1])

    print(name+' simulation at z = ','%.2f' % sim.properties['z'],' and time = ',sim.properties['time'].in_units("Gyr") )

sim.properties
sim.properties['time'].in_units('Gyr')
sim.loadable_keys()
sim.physical_units()
h = sim.halos()
h1 = h[1]
pynbody.analysis.angmom.faceon(h1)

pynbody.plot.sph.image(h1.g,qty='metals', width="500 kpc",filename=name+'_xy_rho_wid500.pdf')
plt.show()

pynbody.plot.sph.image(h1.g,qty='rho', width="50 kpc",filename=name+'_xy_rho_wid50.pdf')
plt.show()


# Constants
m_H = 1.6733 * 10**-24 #g
Z_sun = 0.0142 # (Asplund 2009; https://arxiv.org/pdf/0909.0948.pdf) 

# Want to isolate CGM  
CGM_limit = 10 #kpc
#CGM_limit = 15 #kpc
CGM_gas = h1.g[h1.g['r'].in_units('kpc') > CGM_limit]

print('Total halo mass:',np.sum(h1['mass']))
print('Total gas mass:',np.sum(h1.g['mass']))
print('Total stellar mass:',np.sum(h1.s['mass']))
print('Total CGM gas mass:', np.sum(CGM_gas['mass']))
print('Virial radius:',pynbody.analysis.halo.virial_radius(h1))
#print('Total CGM gas mass in metals:',np.sum(CGM_gas['mass']*CGM_gas['metals'])) # 'metals' *IS* metallicity
CGM_gas['ZoverZsun'] = CGM_gas['metals']/Z_sun
hiZ_CGM_gas = CGM_gas[CGM_gas['ZoverZsun'] >= 0.8]
print('Total CGM gas mass w/ Z > 0.8:',np.sum(hiZ_CGM_gas['mass']))
print('Fraction of Total CGM gas mass that has Z > 0.8:',np.sum(hiZ_CGM_gas['mass'])/np.sum(CGM_gas['mass']))
print('Total CGM gas mass w/ Z > 0.8 and > 20 kpc from center:',np.sum(hiZ_CGM_gas['mass'][hiZ_CGM_gas['r'].in_units('kpc') > 20]))
print('Fraction of CGM gas mass w/ Z > 0.8 that is also > 20 kpc from center:',np.sum(hiZ_CGM_gas['mass'][hiZ_CGM_gas['r'].in_units('kpc') > 20])/np.sum(hiZ_CGM_gas['mass']))

# Plotting Phase Diagrams
x = np.log10(CGM_gas['rho'].in_units('g cm**-3')/m_H)
y = np.log10(CGM_gas['temp'])

# Mass
fig = plt.figure(figsize=(7, 5))
z = CGM_gas['mass']
plt.hexbin(x,y,C=z,reduce_C_function=np.sum,cmap=cm.jet,mincnt=1,bins='log',vmin=5.25,vmax=8.75)
#plt.hexbin(x,y,C=z,reduce_C_function=np.sum,cmap=cm.jet,mincnt=1,bins='log',vmin=1.25,vmax=3.75)
plt.ylabel(r'Log$_{10}$ T/'+str(CGM_gas['temp'].units),size=15)
plt.xlabel(r'Log$_{10}$ n$_H$/cm$^{-3}$',size=15)

from matplotlib.patches import Rectangle
currentAxis = plt.gca()
currentAxis.add_patch(Rectangle((-6,5.2), 8, 0.4, facecolor='Black', alpha=0.2,label='Collisionally Ionized Ovi',hatch='/',edgecolor='Black'))
currentAxis.add_patch(Rectangle((-5,4.8), 1, 0.2, facecolor='Black', alpha=0.5,label='Photoionized Ovi',hatch='|',edgecolor='Black'))

plt.colorbar(label=r'M$_{CGM}$/M$_{\odot}$')
plt.text(-5.5,6.5,name, color='black',size=15)
plt.text(0,6.5,'z = '+str('%.2f' % sim.properties['z']),color='black',size=15)
plt.xlim(-6,2)
plt.ylim(3.5,6.8)
plt.legend(ncol=2,loc=8)
plt.savefig(name+'_phasediagram_CGMat10_mass_'+ts+'_grp.pdf')
plt.show()
plt.clf()

# Metallicity
#fig = plt.figure(figsize=(7, 5))
z = CGM_gas['metals']/Z_sun
plt.hexbin(x,y,C=z,vmin=0.01,vmax=1)#mincnt=1,bins='log',vmin=1.25,vmax=3.75)
plt.ylabel(r'Log$_{10}$ T ('+str(CGM_gas['temp'].units)+')',size=15)
plt.xlabel(r'Log$_{10}$ n$_H$ (cm$^{-3}$)',size=15)

plt.colorbar(label=(r'Z/Z$_{\odot}$'))
plt.text(-5.5,6.5,name, color='black',size=15)
plt.text(0,6.5,'z = '+str('%.2f' % sim.properties['z']),color='black',size=15)
plt.xlim(-6,2)
plt.ylim(3.5,6.8)
plt.savefig(name+'_phasediagram_CGMat10_metallicity_'+ts+'_grp.pdf')
plt.show()
plt.clf()

pynbody.plot.sph.image(CGM_gas,qty='metals', width="500 kpc",filename=name+'_xy_rho_wid500.pdf')
plt.show()

# Radius
#fig = plt.figure(figsize=(7, 5))
z = CGM_gas['r'].in_units('kpc')/int(CGM_gas['r'].max())
plt.hexbin(x,y,C=z,cmap=cm.plasma,vmin=0.1,vmax=1)
plt.ylabel(r'Log$_{10}$ T ('+str(CGM_gas['temp'].units)+')',size=15)
plt.xlabel(r'Log$_{10}$ n$_H$ (cm$^{-3}$)',size=15)

plt.colorbar(label=(r'R/R$_{vir}$'))
plt.text(-5.5,6.7,name, color='black',size=15)
plt.text(0,6.7,'z = '+str('%.2f' % sim.properties['z']),color='black',size=15)
plt.xlim(-6,2)
plt.ylim(3.5,7)
plt.savefig(name+'_phasediagram_CGMat10_Rkpc_'+ts+'_grp.pdf')
#plt.show()
plt.clf()
