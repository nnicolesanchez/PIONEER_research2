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
        sim = pynbody.load('/nobackupp2/nnsanche/NO_BHs/pioneer50h243.1536gst1bwK1_3456/pioneer50h243.1536gst1bwK1.00'+ts)
    elif (str(sys.argv[1]) == 'GM1noBH'):
        sim = pynbody.load('/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM1.1536gst1bwK1_3456/pioneer50h243GM1.1536gst1bwK1.00'+ts)
    elif (str(sys.argv[1]) == 'GM7noBH'):
        sim = pynbody.load('/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM7.1536gst1bwK1_3456/pioneer50h243GM7.1536gst1bwK1.00'+ts) 
    elif (str(sys.argv[1]) == 'GM4noBH'):
        sim = pynbody.load('/nobackupp2/nnsanche/NO_BHs/pioneer50h243GM4.1536gst1bwK1_3456/pioneer50h243GM4.1536gst1bwK1.00'+ts)    

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

# Constants
m_H = 1.6733 * 10**-24 #g
Z_sun = 0.0142 # (Asplund 2009; https://arxiv.org/pdf/0909.0948.pdf) 

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

print('Total CGM gas mass:', np.sum(CGM_gas['mass']))
#print('Total CGM gas mass in metals:',np.sum(CGM_gas['mass']*CGM_gas['metals'])) # 'metals' *IS* metallicity
CGM_gas['ZoverZsun'] = CGM_gas['metals']/Z_sun
hiZ_CGM_gas = CGM_gas[CGM_gas['ZoverZsun'] >= 0.8]
print('Total CGM gas mass w/ Z > 0.8:',np.sum(hiZ_CGM_gas['mass']))
print('Fraction of Total CGM gas mass that has Z > 0.8:',np.sum(hiZ_CGM_gas['mass'])/np.sum(CGM_gas['mass']))
print('Total CGM gas mass w/ Z > 0.8 and > 20 kpc from center:',np.sum(hiZ_CGM_gas['mass'][hiZ_CGM_gas['r'] > 20]))
print('Fraction of CGM gas mass w/ Z > 0.8 that is also > 20 kpc from center:',np.sum(hiZ_CGM_gas['mass'][hiZ_CGM_gas['r'] > 20])/np.sum(hiZ_CGM_gas['mass']))

#print('Total mass Oxygen in CGM:', np.sum(CGM_gas['OxMassFrac']*CGM_gas['mass']),CGM_gas['mass'].units)
#print('Total mass in OVI in CGM:', np.sum(CGM_gas['OxMassFrac']*CGM_gas['mass']*CGM_gas['ovi']))
#print('OVI fractions:',np.average(CGM_gas['ovi']))

# Plotting Phase Diagrams
x = np.log10(CGM_gas['rho'].in_units('g cm**-3')/m_H)
y = np.log10(CGM_gas['temp'])

# Mass
fig = plt.figure(figsize=(7, 5))
z = np.log10(CGM_gas['mass'].in_units('Msol'))
plt.hexbin(x,y,C=z,reduce_C_function=np.sum,cmap=cm.jet,mincnt=1,bins='log',vmin=1.25,vmax=3.75)
plt.ylabel(r'Log$_{10}$ T ('+str(CGM_gas['temp'].units)+')',size=15)
plt.xlabel(r'Log$_{10}$ n$_H$ (cm$^{-3}$)',size=15)

plt.colorbar(label=r'M$_{CGM}$/M$_{\odot}$')
plt.text(-5.5,6.7,name, color='black',size=15)
plt.text(0,6.7,'z = '+str('%.2f' % sim.properties['z']),color='black',size=15)
plt.xlim(-6,2)
plt.ylim(3.5,7)
plt.savefig(name+'_phasediagram_mass_'+ts+'.pdf')
#plt.show()
plt.clf()

# Metallicity
#fig = plt.figure(figsize=(7, 5))
z = CGM_gas['metals']/Z_sun
plt.hexbin(x,y,C=z,vmin=0.01,vmax=1)#mincnt=1,bins='log',vmin=1.25,vmax=3.75)
plt.ylabel(r'Log$_{10}$ T ('+str(CGM_gas['temp'].units)+')',size=15)
plt.xlabel(r'Log$_{10}$ n$_H$ (cm$^{-3}$)',size=15)

plt.colorbar(label=(r'Z/Z$_{\odot}$'))
plt.text(-5.5,6.7,name, color='black',size=15)
plt.text(0,6.7,'z = '+str('%.2f' % sim.properties['z']),color='black',size=15)
plt.xlim(-6,2)
plt.ylim(3.5,7)
plt.savefig(name+'_phasediagram_metallicity_'+ts+'.pdf')
#plt.show()
plt.clf()

# Radius
#fig = plt.figure(figsize=(7, 5))
z = CGM_gas['r'].in_units('kpc')/int(CGM_gas['r'].max())
plt.hexbin(x,y,C=z,cmap=cm.plasma,vmin=0.1,vmax=1)
plt.ylabel(r'Log$_{10}$ T ('+str(CGM_gas['temp'].units)+')',size=15)
plt.xlabel(r'Log$_{10}$ n$_H$ (cm$^{-3}$)',size=15)

plt.colorbar(label=(r'R [kpc]/R$_{vir}$'))
plt.text(-5.5,6.7,name, color='black',size=15)
plt.text(0,6.7,'z = '+str('%.2f' % sim.properties['z']),color='black',size=15)
plt.xlim(-6,2)
plt.ylim(3.5,7)
plt.savefig(name+'_phasediagram_Rkpc_'+ts+'.pdf')
plt.show()
plt.clf()
