import numpy as np
import pynbody

m_h = 1.6733 * 10**-24 # g
galaxies = np.arange(0,300,1)

print('Loading simulation: ROMULUS25')
sim = '/nobackupp2/mtremmel/Romulus/cosmo25/cosmo25p.768sg1bwK1BHe75.006912'
print('Simulation loaded')
s = pynbody.load(sim)
print('Loading halos.')
h = s.halos(dosort=True)
print('Halos loaded.')
#s.physical_units()

M_blue = [9,10,14,15,16,17,18,19,20,21,22,26,27,28,29,30,31,32,33,34,35,36,37,38,42,43,44,45,46,47,48,49,50,51,58,61,63,68,69,70,72,77,79,80,87,89] #no 83,86
M_red  = [11,23,24,41,52,66]
halo   = M_blue + M_red

halomass_array    = []
stellarmass_array = []
id_array          = []
#R_vir_array       = []
#sSFR_array        = []
#O_H_nden          = []
errored = []

# Within 10 kpc - 100 kpc "inner CGM"
total_Mgas        = []
total_Mgas_metals = []

print('Beginning loop.')
for i in range(33,len(halo)):
#    try:
    print('Loading halo:',halo[i])
    h1 = h.load_copy(halo[i])
    h1.physical_units()
    
    pynbody.analysis.angmom.faceon(h1)
    H1_mass      = h1['mass'].sum()
    stellar_mass = h1.s['mass'].sum()
    gas_mass     = h1.g['mass'].sum()
    #R_vir        = pynbody.analysis.halo.virial_radius(h1)    
    
    map_10_100kpc     = ((h1.g['r'].in_units('kpc') > 10) & (h1.g['r'].in_units('kpc') < 100))
    Mgas_10_100kpc    = np.sum(h1.g['mass'][map_10_100kpc])
    Mmetals_10_100kpc = np.sum(h1.g['mass'][map_10_100kpc]*h1.g['metals'][map_10_100kpc])
        
#    except:
#errored.append(i)
#np.savetxt('errored_halo_ids.txt',errored)
#continue

    print('Halo ',galaxies[i],'has mass =',H1_mass,' and stellar mass:',stellar_mass)
    print('Between 10-100 kpc, the total gas mass is:',Mgas_10_100kpc,' and total metal mass is:',Mmetals_10_100kpc)
    #print('Rvir = ',R_vir)#,'and R_max = ',R_max)

    halomass_array.append(H1_mass)
    stellarmass_array.append(stellar_mass)
    id_array.append(i)
#    R_vir_array.append(R_vir)
    total_Mgas.append(Mgas_10_100kpc)
    total_Mgas_metals.append(Mmetals_10_100kpc)

    np.savetxt('COSstellarmass_range_data/ROM_halomasses_z017.txt',halomass_array)
    np.savetxt('COSstellarmass_range_data/ROM_stellarmass_z017.txt',stellarmass_array)
    np.savetxt('COSstellarmass_range_data/ROM_id_z017.txt',id_array)
#    np.savetxt('COSstellarmass_range_data/ROM_Rvir_z017.txt',R_vir_array)
    np.savetxt('COSstellarmass_range_data/ROM_Mtotgas_10_100kpc_z017.txt',total_Mgas)
    np.savetxt('COSstellarmass_range_data/ROM_Mtotmetals_10_100kpc_z017.txt',total_Mgas_metals)
