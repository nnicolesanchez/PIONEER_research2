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
#    halo = [2,4,6,8,10,11,14,16,18,20,22]
#    labels = ['2','4','6','8','10','12','14','16','18','20','22']

#    DD      = [44,32,38,42,47,56,48,27,34,52,31,46]  # Disk Dominated    
#    lab     = [2,4,6,8,10,11,14,16,18,20,22,26,27,29,30,31,32,33,34,35,36,37,38,41,42,43,44,45,46,47,48,49,50,51,52,56]
#    MH_mass = [1.66E+13,1.36E+13,8.6E+12,7.88E+12,5.72E+12,5.12E+12,3.89E+12,3.16E+12,2.66E+12,2.39E+12,2.42E+12,1.54E+12,1.93E+12,1.66E+12,1.50E+12,1.22E+12,1.34E+12,1.38E+12,1.32E+12,1.21E+12,1.09E+12,1.20E+12,9.78E+11,8.05E+11,9.15E+11,7.91E+11,8.42E+11,9.13E+11,8.83E+11,7.63E+11,7.90E+11,7.10E+11,7.21E+11,7.37E+11,7.10E+11,7.59E+11] 


# This array below is for the random sampling version (basically just to make it smaller and fit on one plot)
# See above for 1/2 of high mass galaxies, and ALL MW
# See average_Ofrac... script for the 50low mass galaxies I averaged
    lab     = [2,4,6,8,10,11,14,16,18,20,22,26,30,32,34,36,38,42,44,46,48,50,52,56,60,70,80,90,100,110,120,135,145,155,170,180,190,200,210,220,230,240,250,260,270,280,300]
    MH_mass = [1.66E+13,1.36E+13,8.6E+12,7.88E+12,5.72E+12,5.12E+12,3.89E+12,3.16E+12,2.66E+12,2.39E+12,2.42E+12,1.54E+12,1.50E+12,1.34E+12,1.32E+12,1.09E+12,9.78E+11,9.15E+11,8.42E+11,8.83E+11,7.90E+11,7.21E+11,7.10E+11,7.59E+11,6.432556e+11,4.111969e+11,3.477979e+11,3.461157e+11,1.9526505e+11,2.4981697e+11,2.2209713e+11,1.5930106e+11,1.9665586e+11,1.6736022e+11,1.6696692e+11,1.9500245e+11,2.0052886e+11,1.432943e+11,1.4283268e+11,1.4880883e+11,1.3410762e+11,1.0599081e+11,1.3814106e+11,1.2795174e+11,5.0811847e+10,1.01250515e+11,1.1392372e+11] # MH mass


    MH_mass_sort = np.argsort(np.array(MH_mass))
    print(lab)
    lab = np.array(lab)
    lab     = lab[MH_mass_sort]
    MH_mass = np.sort(MH_mass)
    
#    lab     = lab[::-1]
#    print(lab)
#    MH_mass = MH_mass[::-1]

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

text_place = -0.7
for m in range(len(MH_mass)):
#    print(lab[m],"%.1f" % np.log10(MH_mass[m]),MH_mass[m])
    plt.text(text_place,1.1,str("%.1f" % np.log10(MH_mass[m])),color='Black',rotation=90)
    text_place = text_place + 2

plt.ylim(0,1.4)
plt.legend(ncol=4)
plt.savefig('OPPEN_fig10_ROMgxys_plushighmass.pdf')
plt.show()


