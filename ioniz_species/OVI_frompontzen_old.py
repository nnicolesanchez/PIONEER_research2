import matplotlib.pyplot as plt
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
    print(f.gas['rho'].in_units('g cm**-3')*ovi*f.gas['OxMassFrac']/(16*m_p))
    return f.gas['rho'].in_units('g cm**-3')*ovi*f.gas['OxMassFrac']/(16*m_p)

    
k = 4
sim = ['/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.00']
labels = ['P0','GM1','GM4','GM5','GM6','GM7']
ts = np.loadtxt('../'+labels[k]+'/timesteps.txt',dtype=str)
t = len(ts)-1

f = pynbody.load(sim[k]+ts[t])

pynbody.analysis.halo.center(f.star)

f.physical_units()
h = f.halos()
h1 = h[1]
#print(h1.properties['Rvir'])
#print(f.gas['rho'].units)
#print(ovi)
pynbody.analysis.angmom.faceon(h1)


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


# From Pontzen's code
pynbody.plot.sph.image(CGM_gas,qty='N_OVI',width='500 kpc',cmap='magma',show_cbar=False,units="g cm^-2",vmin=1e13,vmax=1e15)

plt.text(200,200,str(labels[k]),color='White')
plt.colorbar(label=r'N$_{OVI}$ cm$^{-2}$')
plt.savefig(str(labels[k])+'_xvsy_Novi.pdf')   
plt.show()





quit()
############################
# CALCULATE COLUMN DENSITY #
############################
z_max = np.max(CGM_gas['z'])
z_min = np.min(CGM_gas['z'])
L = z_max - z_min 
print(L,L.units)
print(z_min,z_max)

N_gas = CGM_gas['rho']*L
print(N_gas.units)

ovi = pynbody.analysis.ionfrac.calculate(CGM_gas,ion='ovi',mode='new') 
m_p = 1.6726 * 10**-24 #g
#print(ovi)
#print(f.gas['OxMassFrac'])
N_OVI = N_gas.in_units('g cm**-2')*ovi*CGM_gas['OxMassFrac']/(16*m_p)
print('Total mass Oxygen in CGM:', np.sum(CGM_gas['OxMassFrac']*CGM_gas['mass']),CGM_gas['mass'].units)
print('Total mass in OVI in CGM:', np.sum(CGM_gas['OxMassFrac']*CGM_gas['mass']*ovi))
print('OVI Column Density: ',N_OVI,N_OVI.units)

r = (CGM_gas['x']**2 + CGM_gas['y']**2)**(0.5)
print(r,CGM_gas['x'].units)
print(np.min(r))

n = np.arange(0,240,10)
print(n)
N_OVI_avg = []
N_OVI_median = []
N_OVI_hist = []
sum = 0 
for i in range(len(n)-1):
    CGM_gas_mask = [r >= n[i]] and [r < n[i+1]]
    N_OVI_bin = N_OVI[CGM_gas_mask]
    N_OVI_avg.append(np.mean(N_OVI_bin))
    N_OVI_median.append(np.median(N_OVI_bin))
    N_OVI_hist.append(np.sum(N_OVI_bin))
    
n_mid = n + 5
print(n_mid)
plt.plot(n_mid[:-1],np.log10(N_OVI_avg),label='Mean')
#plt.ylabel(r'N$_{OVI}$',label='Mean')
#plt.xlabel(r'$b$ [kpc]')
#plt.show()
plt.plot(n_mid[:-1],np.log10(N_OVI_median),label='Median')
plt.ylabel(r'N$_{OVI}$')
plt.xlabel(r'$b$ [kpc]')
plt.ylim(14,17.5)
plt.xlim(0,200)
plt.text(140,17.25,str(labels[k]))
plt.legend()
plt.savefig(labels[k]+'_NOVI_b.pdf')
plt.show()


#quit()
##########################################
# PLOT COLUMN DENSITY AS FUNCTION OF X,Y #
##########################################
# From Pontzen's code
pynbody.plot.sph.image(CGM_gas,qty='N_OVI',width='500 kpc',cmap='magma',show_cbar=False,units="g cm^-2",vmin=1e13,vmax=1e15)

plt.text(200,200,str(labels[i]),color='White')
plt.colorbar(label=r'N$_{OVI}$ cm$^{-2}$')
#plt.savefig(str(labels[i])+'_OVI.pdf')
plt.show()












