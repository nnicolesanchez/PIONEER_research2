# This script recreates the plot from Oppenheimer 2016 Figure 10
# Oxygen mass fractions as a function of halo mass (bar plot)
# Using hdf5 table with ion fractions from Trident


# N. Nicole Sanchez -- Created: June 27, 2018
# Univ. of Wash, S  -- Edited: June 27, 2018
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib
import sys

if len(sys.argv) == 1:
    print('No selection made. Options are "BH" or "noBH"')
    print('Syntax: Ofrac_barplot.py BH')
    quit()
elif (str(sys.argv[1]) == 'BH'):
    print('Loading in ROM')
    DD      = [44,32,38,42,47,56,48,27,34,52,31,46]  # Disk Dominated    
    lab     = [2,4,6,8,10,12,14,16,18,20,22,26,27,29,30,31,32,33,34,35,36,37,38,41,42,43,44,45,46,47,48,49,50,51,52,56]
    lab     = lab[::-1]
    print(lab)
    #halo = [2,4,6,8,10,11,14,16,18,20,22]
#    labels = ['2','4','6','8','10','12','14','16','18','20','22']

    MH_mass = [1.66E+13,1.36E+13,8.6E+12,7.88E+12,5.72E+12,5.12E+12,3.89E+12,3.16E+12,2.66E+12,2.39E+12,2.42E+12,1.54E+12,1.93E+12,1.66E+12,1.50E+12,1.22E+12,1.34E+12,1.38E+12,1.32E+12,1.21E+12,1.09E+12,1.20E+12,9.78E+11,8.05E+11,9.15E+11,7.91E+11,8.42E+11,9.13E+11,8.83E+11,7.63E+11,7.90E+11,7.10E+11,7.21E+11,7.37E+11,7.10E+11,7.59E+11] # log of MH mass
    MH_mass = MH_mass[::-1]

#    print('Loading in GM runs with BH physics.')
#    lab     = ['P0','GM1','GM7','GM4']
#    NEW_lab = ['P0','GM1','GM2','GM3']
#elif (str(sys.argv[1]) == 'noBH'):
#    print('Loading in GM runs with noBH physics.')
#    lab     = ['P0noBH','GM1noBH','GM7noBH','GM4noBH']
#    NEW_lab = ['P0noBH','GM1noBH','GM2noBH','GM3noBH']
else: 
    print('No valid option input. Goodbye.')
    quit()

ion_labels = ['oi','oii','oiii','oiv','ov','ovi','ovii','oviii']

# Attempting to match color scheme of Opp 2016 Fig 10
colors = ['darkmagenta','magenta','red','orange','gold','lime','cyan','blue']

# Setup bar plot parameters
N   = len(lab)   # Number of bars
indec = np.arange(0,2*N,2)
wid = 1.5

ion_array  = []
# Read in the OxMass and Ion Fractions to calculate the amount of 
# mass in every ion and total oxygen mass

plt.figure(num=None, figsize=(10, 4))

for i in range(len(ion_labels)):
    print('Making ion array for:',ion_labels[i])
    ion_mass_frac = []

    for s in range(len(lab)):
        print(s)
        Omass    = np.loadtxt('/home1/nnsanche/ROMULUS_research/ALLcolumndensities/Omass_z0_17/'+str(lab[s])+'_Omass_3456.np')
        ion_frac = np.loadtxt('/home1/nnsanche/ROMULUS_research/ALLcolumndensities/ion_frac_z0_17/'+str(lab[s])+'_'+ion_labels[i]+'frac_3456.np') 

#        Omass      = np.loadtxt('../profiles_plot7/'+lab[s]+'/'+lab[s]+'_Omass_3456.np')
#        ion_frac   = np.loadtxt('../profiles_plot7/'+lab[s]+'/'+lab[s]+'_'+ion_labels[i]+'frac_3456.np')
        ion_mass_frac.append(np.nansum(ion_frac*Omass)/np.nansum(Omass))
        print(ion_mass_frac)
        
    if i == 0 :
        print('Plotting oi')
        plt.bar(indec, ion_mass_frac,width=wid,label=ion_labels[i],color=colors[i],edgecolor='Black')
        ion_butt = np.array(ion_mass_frac)
    else :
        print('Plotting other ions')
        plt.bar(indec, ion_mass_frac,width=wid,label=ion_labels[i],color=colors[i],edgecolor='Black',bottom=ion_butt)
        ion_butt = np.array(ion_butt) + np.array(ion_mass_frac)

plt.ylabel(r'Oxygen Fraction (CGM <R$_{200}$)')

text_place = -0.6
for m in range(len(MH_mass)):
#    print(lab[m],"%.1f" % np.log10(MH_mass[m]),MH_mass[m])
    plt.text(text_place,0.1,str("%.1f" % np.log10(MH_mass[m])),color='White',rotation=90)
    text_place = text_place + 2

plt.legend(ncol=4)
plt.savefig('OPPEN_fig10_ROMgxys_plushighmass.pdf')
plt.show()


