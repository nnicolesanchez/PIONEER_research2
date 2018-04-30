# This scripts creates a plot comparable to that of 
# Figure 8 in Tumlinson, Peebles, Werk 2017 review.

# N. Nicole Sanchez -- Pleiades:~/sanchenn/PIONEER*/gastempfractions/T_vs_Mstar.py
# UW, Seattle       -- Created: February 26, 2018
from pynbody.analysis import profile
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import numpy as np
import os.path
import pynbody


#k      = 0
labels = ['P0','GM1','GM4','GM5','GM6','GM7']
colors = ['DodgerBlue','SteelBlue','FireBrick','IndianRed','Salmon','Orange']
sims   = ['/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM1.1536gs1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.00','/nobackupp8/ambrook2/fgoverna_pleiades_p8_files/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.00','/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.00']

gas_labels = ['stars','ISM','t < 10^4','low ionization','OVI-traced','X-ray-traced']
gas_colors = ['darkred','skyblue','salmon','purple','yellowgreen','gold']

for k in range(len(sims)):
    ts = np.loadtxt('../'+labels[k]+'/timesteps.txt',dtype=str)
    t  = len(ts)-1
    if k == 5:
        t  = len(ts)-2 # Since GM7 still has a problem with timestep 4096

    # Load in Sim
    sim = pynbody.load(sims[k]+ts[t])
    print('Loading sim:',sims[k],' at timestep:',ts[t])
    
    # Focus on Main Halo
    h  = sim.halos()
    h1 = h[1]
    pynbody.analysis.halo.center(h1,mode='com')
    pynbody.analysis.angmom.faceon(h1)
    sim.physical_units()
        
    # Total Mass of the Halo
    TOT_halo_mass = h1.g['mass'].sum()+h1.s['mass'].sum()
    
    # Total Mass of Stars
    TOT_stellar_mass = h1.s['mass'].sum()
    
    ###################
    # ISOLATE CGM GAS #   
    ###################
    # Isolate and remove disk stars within radius 0-10 kpc & vertically 10 kpc 
    r_max = 10  # kpc
    twenty_kpc_incm = 6.171*(10**22)
    
    Rg_d = ((h1.g['x'].in_units('kpc'))**2. + (h1.g['y'].in_units('kpc'))**2. + (h1.g['z'].in_units('kpc'))**2.)**(0.5)
    disk_gas_xyzmax =  (Rg_d < r_max)
    disk_gas_mask = disk_gas_xyzmax #& disk_gas_zmax
    ISM_gas = h1.g[disk_gas_mask] #& disk_gas_zmax]
    
    # Total ISM gas
    TOT_ISM_gas = ISM_gas['mass'].sum()
    
    #########################
    # Split CGM into Phases #
    #########################
    # (useful?) T < 10^4
    # Cold Gas: 10^4 < T < 10^5
    # Warm Gas: 10^5 < T < 10^6
    # Hot Gas : T > 10^6
    CGM_gas  = h1.g[~disk_gas_mask]
        
    COLDER_CGM_gas = CGM_gas[(CGM_gas['temp'] <= 10**4)]
    COLD_CGM_gas = CGM_gas[(CGM_gas['temp'] > 10**4) & (CGM_gas['temp'] <= 10**5)]
    WARM_CGM_gas = CGM_gas[(CGM_gas['temp'] > 10**5) & (CGM_gas['temp'] <= 10**6)] 
    HOT_CGM_gas  = CGM_gas[(CGM_gas['temp'] > 10**6)]
    
    TOT_colderCGMgas = COLDER_CGM_gas['mass'].sum()
    TOT_coldCGMgas = COLD_CGM_gas['mass'].sum()
    TOT_warmCGMgas = WARM_CGM_gas['mass'].sum()
    TOT_hotCGMgas  = HOT_CGM_gas['mass'].sum()

    #totcold    = (h258galcold, h258bhcold)
    #totclumpy  = (h258galclumpy, h258bhclumpy)
    #totshocked = (h258galshocked, h258bhshocked)
    #totcoldandclumpy = np.array(totclumpy) + np.array(totcold)
    #N   = 2
    #ind = np.arange(N)
    #width = 0.8
    #cold    = plt.bar(ind, totcold, width=width,color='DodgerBlue',edgecolor=outline, bottom=totclumpy)
    #clumpy  = plt.bar(ind, totclumpy, width=width,color='LimeGreen',edgecolor=outline)
    #shocked = plt.bar(ind, totshocked, width=width,color='Red',edgecolor=outline, bottom=totcoldandclumpy)
    
    stellar_frac = TOT_stellar_mass/TOT_halo_mass
    ISM_frac = TOT_ISM_gas/TOT_halo_mass
    CGMcolder_frac = TOT_colderCGMgas/TOT_halo_mass
    CGMcold_frac = TOT_coldCGMgas/TOT_halo_mass
    CGMwarm_frac = TOT_warmCGMgas/TOT_halo_mass
    CGMhot_frac = TOT_hotCGMgas/TOT_halo_mass
    print(stellar_frac,ISM_frac,CGMcolder_frac,CGMcold_frac,CGMwarm_frac,CGMhot_frac)
    
    wid = 0.03

    if k == 0:
        plt.bar(np.log10(TOT_stellar_mass),TOT_stellar_mass/TOT_halo_mass,color=colors[k],width=wid,label=gas_labels[0])#label=labels[k])
        plt.bar(np.log10(TOT_stellar_mass),TOT_ISM_gas/TOT_halo_mass,color=gas_colors[1],bottom=TOT_stellar_mass/TOT_halo_mass,width=wid,label=gas_labels[1])
        plt.bar(np.log10(TOT_stellar_mass),TOT_colderCGMgas/TOT_halo_mass,color=gas_colors[2],bottom=(TOT_stellar_mass+TOT_ISM_gas)/TOT_halo_mass,width=wid,label=gas_labels[2])
        plt.bar(np.log10(TOT_stellar_mass),TOT_coldCGMgas/TOT_halo_mass,color=gas_colors[3],bottom=(TOT_stellar_mass+TOT_ISM_gas+TOT_colderCGMgas)/TOT_halo_mass,width=wid,label=gas_labels[3])
        plt.bar(np.log10(TOT_stellar_mass),TOT_warmCGMgas/TOT_halo_mass,color=gas_colors[4],bottom=(TOT_stellar_mass+TOT_ISM_gas+TOT_coldCGMgas+TOT_colderCGMgas)/TOT_halo_mass,width=wid,label=gas_labels[4])
        plt.bar(np.log10(TOT_stellar_mass),TOT_hotCGMgas/TOT_halo_mass,color=gas_colors[5],bottom=(TOT_stellar_mass+TOT_ISM_gas+TOT_coldCGMgas+TOT_colderCGMgas+TOT_warmCGMgas)/TOT_halo_mass,width=wid,label=gas_labels[5])
        plt.bar(np.log10(TOT_stellar_mass),TOT_stellar_mass/TOT_halo_mass,color=colors[k],width=wid,label=labels[k])
    else:
        plt.bar(np.log10(TOT_stellar_mass),TOT_stellar_mass/TOT_halo_mass,color=colors[k],width=wid,label=labels[k])
        plt.bar(np.log10(TOT_stellar_mass),TOT_ISM_gas/TOT_halo_mass,color=gas_colors[1],bottom=TOT_stellar_mass/TOT_halo_mass,width=wid)
        plt.bar(np.log10(TOT_stellar_mass),TOT_colderCGMgas/TOT_halo_mass,color=gas_colors[2],bottom=(TOT_stellar_mass+TOT_ISM_gas)/TOT_halo_mass,width=wid)
        plt.bar(np.log10(TOT_stellar_mass),TOT_coldCGMgas/TOT_halo_mass,color=gas_colors[3],bottom=(TOT_stellar_mass+TOT_ISM_gas+TOT_colderCGMgas)/TOT_halo_mass,width=wid)
        plt.bar(np.log10(TOT_stellar_mass),TOT_warmCGMgas/TOT_halo_mass,color=gas_colors[4],bottom=(TOT_stellar_mass+TOT_ISM_gas+TOT_coldCGMgas+TOT_colderCGMgas)/TOT_halo_mass,width=wid)
        plt.bar(np.log10(TOT_stellar_mass),TOT_hotCGMgas/TOT_halo_mass,color=gas_colors[5],bottom=(TOT_stellar_mass+TOT_ISM_gas+TOT_coldCGMgas+TOT_colderCGMgas+TOT_warmCGMgas)/TOT_halo_mass,width=wid)


plt.ylabel(r'Cumulative Fraction, M/M$_{h}$')
plt.xlabel(r'log M$_{*}$/M$_{\odot}$')
plt.ylim(0,1.2)
plt.xlim(9.8,11)
plt.legend(ncol=2)
plt.savefig('CumulativeGasFractionsvsStellarMass.pdf')
plt.show()
    
