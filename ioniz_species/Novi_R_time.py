import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import os.path
import seaborn as sns
import pandas as pd
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

# Just using k = 1 and k = 2, for GM1 & GM4 for now
k = 2
## MOVED FILES FROM FABIO TO ALYSON BROOKS: /nobackupp8/ambrook2/fgoverna_pleiades_p8_files
sim = ['/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.00','/nobackupp8/fgoverna/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.004096']
labels = ['P0','GM1','GM4','GM5','GM6','GM7']
colors = sns.cubehelix_palette(8)
Z_sun = 0.0142 # (Asplund 2009; https://arxiv.org/pdf/0909.0948.pdf)
print('LOADING SIM:',labels[k])

ts = np.loadtxt('../'+labels[k]+'/timesteps.txt',dtype=str)
for t in range(len(ts)):
#    t = len(ts)-1
    print('Loading sim:',sim[k],' at timestep:',ts[t])
    if os.path.isfile(labels[k]+'_ionizb_'+ts[t]+'_xdata.np'): 
        continue
        print('This plot has already been made. Skipping this snapshot.x')

    ######################
    # READ IN SIMULATION #
    ######################
    f = pynbody.load(sim[k]+ts[t])
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

    CGM_temp = np.array(CGM_gas['temp'])
    transit_temp = 10**5.2
    CGM_gas = CGM_gas[CGM_temp < transit_temp]
    print('LOOKING AT CGM WITH TEMP < 10^5.2')

    #########################
    # CALCULATE OVI DENSITY #
    #########################
    # Ionization fraction of OVI compare to total Oxygen
    ovi = pynbody.analysis.ionfrac.calculate(CGM_gas,ion='ovi',mode='new') 
    m_p = 1.6726 * 10**-24 #g
    
    # OVI density = total CGM gas density * fraction of oxygen * fraction of OVI / mass of oxygen
    OVI = CGM_gas['rho'].in_units('g cm**-3')*ovi*CGM_gas['OxMassFrac']/(16*m_p)
    print('Total mass in CGM:', np.sum(CGM_gas['mass']))
    if np.sum(CGM_gas['mass']) == 0.0:
        continue
    print('Total mass in metals:',np.sum(CGM_gas['mass']*CGM_gas['metals'])) # 'metals' *IS* metallicity
    print('Total mass Oxygen in CGM:', np.sum(CGM_gas['OxMassFrac']*CGM_gas['mass']),CGM_gas['mass'].units)
    print('Total mass in OVI in CGM:', np.sum(CGM_gas['OxMassFrac']*CGM_gas['mass']*ovi))
    #print('OVI Density: ',OVI,OVI.units)
    print('OVI fractions:',np.average(ovi))
    


    #############################################################
    # DIVIDE PARTICLES INTO SHELLS & CALCULATE COLUMN DENSITIES #
    #############################################################
    # In shells of 10 kpc, take an average density and line of sight L for each shell
    shell_bounds = np.arange(0,275,10)
    CGM_r = (CGM_gas['x']**2 + CGM_gas['y']**2)**(0.5)
    R_vir = int(np.max(CGM_r))
    print('R_vir',int(np.max(CGM_r)),CGM_gas['x'].units)
    

    #fig = plt.figure(figsize=(8, 8))
    #ax1 = fig.add_subplot(221)
    #ax2 = fig.add_subplot(222)
    CGM_Novi = []
    pathlength = []
    for i in range(len(shell_bounds)-1):
        shell = CGM_gas[(np.abs(CGM_r) > shell_bounds[i]) & (np.abs(CGM_r) < shell_bounds[i+1])]
        shell_ovi_frac = ovi[(np.abs(CGM_r) > shell_bounds[i]) & (np.abs(CGM_r) < shell_bounds[i+1])]
        shell_OVI_rho = shell['rho'].in_units('g cm**-3')*shell_ovi_frac*shell['OxMassFrac']/(16*m_p)
        #    pynbody.plot.sph.image(shell,qty="rho",units="g cm^-3",width=300,z_camera=)
        
        # Test some stuff to make sure I'm doing my path length measurements right
        shell_x = np.array(shell['x'])
        shell_y = np.array(shell['y'])
        shell_z = np.array(shell['z'])
        
        if len(shell_x) == 0:
            break

        #ax1.plot(shell_x,shell_y,'.')#color=sns.cubehelix_palette(8)[i])
        j = 1 - i/26.
        #print(j)
        #ax2.plot(shell_x,shell_z,'.',alpha=j)
        #print(np.min(shell_z),np.max(shell_z))
        
        avg_shell_OVI_rho = np.average(shell_OVI_rho)
        #print(avg_shell_OVI_rho,shell_OVI_rho.units)
        
        shell_z = np.max(shell['z'].in_units('cm')) - np.min(shell['z'].in_units('cm'))
    
        if i <= 1 :
            CGM_Novi.append(avg_shell_OVI_rho*(shell_z - twenty_kpc_incm))
            print('Remove 20 kpc in z because of empty center within 10 kpc radius')
        else:
            CGM_Novi.append(avg_shell_OVI_rho*shell_z) # Using shell_z to underestimate gas
            pathlength.append(shell_z/(3.086*10**21))

    #plt.show()
    #print(CGM_Novi)

    b_impact = shell_bounds + 5
    #print(b_impact)

    pathlength_kpc = pathlength #kpc

    #plt.plot(shell_bounds[:-3],pathlength_kpc)
    #plt.show()
    #quit()

    #print(float(ts[t]))
    #if (float(ts[t]) >= 2816.):
        #COS = pd.read_csv('COShalo_obs.txt',header=0,delim_whitespace=True,comment='#',index_col=False)
        #COS_ID = COS['ID']
        #COS_OVI = COS['LogN0VI']
        #COS_Rkpc = COS['Rho(kpc)'][COS_OVI != 0]
        #COS_OVI = COS_OVI[COS_OVI != 0]
        #print(COS)
        #plt.plot(COS_Rkpc[COS['RED?'] == 'yes'],COS_OVI[COS['RED?'] == 'yes'],marker='.',color='Red',linestyle=' ',label='COS Ellipticals')
        #plt.plot(COS_Rkpc[COS['RED?'] == 'no'],COS_OVI[COS['RED?'] == 'no'],marker='.',color='DodgerBlue',linestyle=' ',label='COS Spirals')
        
        
    np.savetxt(labels[k]+'_ionizb_'+ts[t]+'_xdata.np',shell_bounds[0:len(CGM_Novi)])
    np.savetxt(labels[k]+'_ionizNovi_'+ts[t]+'_ydata.np',CGM_Novi)
#    plt.plot(shell_bounds[0:len(CGM_Novi)],np.log10(CGM_Novi),marker='.',label=labels[k])
#    plt.ylabel(r'N$_{OVI}$ [cm$^{-2}$]')
#    plt.xlabel(r'$r$ [kpc]')
#    plt.ylim(13,16)
#    plt.xlim(-10,270)
#    plt.text(-5,15.75,'Collisional OVI (T > 10^${5.2}$) at '+str(labels[k]))
#    plt.legend()
#    plt.savefig(labels[k]+'_collNOVI_b_'+ts[t]+'.pdf')
#    plt.show()
#    quit()

quit()
##########################################
# PLOT COLUMN DENSITY AS FUNCTION OF X,Y #
##########################################
# From Pontzen's code
pynbody.plot.sph.image(CGM_gas,qty='N_OVI',width='500 kpc',cmap='magma',show_cbar=False,units="g cm^-2",vmin=1e13,vmax=1e15)

plt.text(200,200,str(labels[i]),color='White')
plt.colorbar(label=r'N$_{OVI}$ cm$^{-2}$')
#plt.savefig(str(labels[i])+'_OVI.pdf')
plt.show()













