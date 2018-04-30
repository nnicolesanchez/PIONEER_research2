import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import seaborn as sns
import numpy as np
import pynbody

@pynbody.derived_array
#def rhoOVI(f):
#    ovi = pynbody.analysis.ionfrac.calculate(f,ion='ovi',mode='new')
#    return f.gas['rho']*ovi*f.gas['OxMassFrac']  # Original from Pontzen

def N_OVI(f):
    ovi = pynbody.analysis.ionfrac.calculate(f,ion='ovi',mode='new')
    m_p = 1.6726 * 10**-24 #g
    #print(ovi)
    #print(f.gas['OxMassFrac'])
    #print(f.gas['rho'].in_units('g cm**-3')*ovi*f.gas['OxMassFrac']/(16*m_p))
    return f.gas['rho'].in_units('g cm**-3')*ovi*f.gas['OxMassFrac']/(16*m_p)

    
k = 1
sim = ['/nobackupp8/fgoverna/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.004096','/nobackupp8/fgoverna/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.004096','/nobackupp8/fgoverna/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.004096','/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.004096','/nobackupp8/fgoverna/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.004096','/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.004096']
labels = ['P0','GM1','GM4','GM5','GM6','GM7']
colors = sns.cubehelix_palette(8)
Z_sun = 0.0142 # (Asplund 2009; https://arxiv.org/pdf/0909.0948.pdf)
print('LOADING SIM:',labels[k])

######################
# READ IN SIMULATION #
######################
f = pynbody.load(sim[k])
pynbody.analysis.halo.center(f.star)
f.physical_units()
h = f.halos()
h1 = h[1]
pynbody.analysis.angmom.faceon(h1)

###################
# ISOLATE CGM GAS #   
###################
# Isolate and remove disk stars within radius 0-10 kpc & vertically 10 kpc 
r_max = 10  # kpc
twenty_kpc_incm = 6.171*(10**22)
#z_max = 10 #4 # kpc

Rg_d = ((h1.g['x'].in_units('kpc'))**2. + (h1.g['y'].in_units('kpc'))**2. + (h1.g['z'].in_units('kpc'))**2.)**(0.5)
disk_gas_xyzmax =  (Rg_d < r_max)
#disk_gas_zmax  = (h1.g['z'].in_units('kpc') < z_max) & (h1.g['z'].in_units('kpc') > -z_max)
disk_gas_mask = disk_gas_xyzmax #& disk_gas_zmax
disk_gas = h1.g[disk_gas_mask] #& disk_gas_zmax]
CGM_gas  = h1.g[~disk_gas_mask]


#########################
# CALCULATE OVI DENSITY #
#########################
# Ionization fraction of OVI compare to total Oxygen
ovi = pynbody.analysis.ionfrac.calculate(CGM_gas,ion='ovi',mode='new') 
m_p = 1.6726 * 10**-24 #g

# OVI density = total CGM gas density * fraction of oxygen * fraction of OVI / mass of oxygen
OVI = CGM_gas['rho'].in_units('g cm**-3')*ovi*CGM_gas['OxMassFrac']/(16*m_p)
print('Total mass in CGM:', np.sum(CGM_gas['mass']))
print('Total mass in metals:',np.sum(CGM_gas['mass']*CGM_gas['metals'])) # 'metals' *IS* metallicity
print('Total mass Oxygen in CGM:', np.sum(CGM_gas['OxMassFrac']*CGM_gas['mass']),CGM_gas['mass'].units)
print('Total mass in OVI in CGM:', np.sum(CGM_gas['OxMassFrac']*CGM_gas['mass']*ovi))
#print('OVI Density: ',OVI,OVI.units)
print('OVI fractions:',np.average(ovi))

CGM_temp = np.array(CGM_gas['temp'])
#print(len(CGM_temp[CGM_temp < 10000]))
#print(np.min(CGM_temp))

#plt.plot(np.log10(CGM_temp[10000:20000]),ovi[10000:20000],'.')
#plt.plot(np.log10(CGM_temp[CGM_temp < 10000]),ovi[CGM_temp < 10000],'.')
#plt.ylabel(r'f$_{OVI}$')
#plt.xlabel('log(T [K])')
#plt.savefig(str(labels[k])+'_fracOVI_temp.pdf')
#plt.show()

#######################################################
# Color code f_ovi vs temp by density and metallicity #
#######################################################
# METALLICITY
CGM_metallicity = CGM_gas['metals']
fig = plt.figure(figsize=(7, 5))
plt.hist2d(np.log10(CGM_temp),ovi,(100,100),weights=CGM_gas['metals'],cmap=cm.jet,norm=mpl.colors.LogNorm())
plt.ylabel(r'$f_{OVI}$')
plt.xlabel(r'Log(T [K])')
plt.colorbar(label=(r'Z/Z$_{\odot}$'),ticklabels=True)
plt.text(3.9,0.185,labels[k], color='midnightblue',size=12)
#plt.text(0.7,6.7,'z = 0',color='midnightblue',size=12)
plt.xlim(3.7,6.7)
plt.ylim(0,0.2)
plt.savefig(labels[k]+'_fracovi_temp_metallicity.pdf')
plt.show()

quit()

# DENSITY
fig = plt.figure(figsize=(7, 5))
plt.hist2d(np.log10(CGM_temp),ovi,(100,100),weights=(CGM_gas['rho'].in_units('g cm**-3')),cmap=cm.jet,norm=mpl.colors.LogNorm())
plt.ylabel(r'$f_{OVI}$')
plt.xlabel(r'Log(T [K])')
plt.colorbar(label=str('g cm$^{-3}$'))
plt.text(3.9,0.185,labels[k], color='midnightblue',size=12)
#plt.text(0.7,6.7,'z = 0',color='midnightblue',size=12) 
plt.xlim(3.7,6.7)
plt.ylim(0,0.2)
plt.savefig(labels[k]+'_fracovi_temp_density.pdf')
plt.show()


######################
# DENSITY HISTOGRAMS #
######################
# T < 10^5
CGM_T_lt10_5 = CGM_gas[CGM_temp < 10**5.]
plt.hist(CGM_T_lt10_5)
plt.show()



# f_OVI < 0.05
CGM_fOVI_lt0_05 = CGM_gas[ovi < 0.05]
plt.hist(CGM_fOVI_lt0_05)
plt.show()






















