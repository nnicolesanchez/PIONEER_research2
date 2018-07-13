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
#    lab     = [2,4,6,8,10,11,14,16,18,20,22,26,27,29,30,31,32,33,34,35,36,37,38,41,42,43,44,45,46,47,48,49,50,51,52,56]


    lab     = [2,4,6,8,10,11,14,16,18,20,22,26,27,29,30,31,32,33,34,35,36,37,38,41,42,43,44,45,46,47,48,49,50,51,52,56,60,70,80,90,95,100,105,110,115,120,125,135,140,145,150,155,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240]#,245,250,255,260,265,270,275,280,295,300]
    lab     = lab[::-1]
    print(lab)
    #halo = [2,4,6,8,10,11,14,16,18,20,22]
#    labels = ['2','4','6','8','10','12','14','16','18','20','22']

#    MH_mass = [1.66E+13,1.36E+13,8.6E+12,7.88E+12,5.72E+12,5.12E+12,3.89E+12,3.16E+12,2.66E+12,2.39E+12,2.42E+12,1.54E+12,1.93E+12,1.66E+12,1.50E+12,1.22E+12,1.34E+12,1.38E+12,1.32E+12,1.21E+12,1.09E+12,1.20E+12,9.78E+11,8.05E+11,9.15E+11,7.91E+11,8.42E+11,9.13E+11,8.83E+11,7.63E+11,7.90E+11,7.10E+11,7.21E+11,7.37E+11,7.10E+11,7.59E+11] # log of MH mass
#    MH_mass = MH_mass[::-1]

else: 
    print('No valid option input. Goodbye.')
    quit()

ion_labels = ['oi','oii','oiii','oiv','ov','ovi','ovii','oviii']

# Attempting to match color scheme of Opp 2016 Fig 10
colors = ['darkmagenta','magenta','red','orange','gold','lime','cyan','blue']

# Setup bar plot parameters
N   = 3 #len(lab)   # Number of bars
indec = np.arange(0,2*N,2)
wid = 1.0

ion_array = []
# Read in the OxMass and Ion Fractions to calculate the amount of 
# mass in every ion and total oxygen mass

for i in range(len(ion_labels)):
    lo_ion_frac = []
    mid_ion_frac = []
    hi_ion_frac  = []
    med_array = []
    print('Making ion array for:',ion_labels[i])

    print('Looking at Low-mass')
    for s in range(len(lab)-36):
        print('in galaxies',lab[s])

        Omass    = np.loadtxt('/home1/nnsanche/ROMULUS_research/ALLcolumndensities/Omass_z0_17/'+str(lab[s])+'_Omass_3456.np')
        ion_frac = np.loadtxt('/home1/nnsanche/ROMULUS_research/ALLcolumndensities/ion_frac_z0_17/'+str(lab[s])+'_'+ion_labels[i]+'frac_3456.np')

        lo_ion_frac.append(np.nansum(ion_frac*Omass)/np.nansum(Omass))
        print('ion fraction =',np.nansum(ion_frac*Omass)/np.nansum(Omass))

    med_array.append(np.median(lo_ion_frac))

    print('Now looking at MW-mass')
    for s in range(len(lab)-36,len(lab)-11):
        print('in galaxies',lab[s])

        Omass    = np.loadtxt('/home1/nnsanche/ROMULUS_research/ALLcolumndensities/Omass_z0_17/'+str(lab[s])+'_Omass_3456.np')
        ion_frac = np.loadtxt('/home1/nnsanche/ROMULUS_research/ALLcolumndensities/ion_frac_z0_17/'+str(lab[s])+'_'+ion_labels[i]+'frac_3456.np') 
        
        mid_ion_frac.append(np.nansum(ion_frac*Omass)/np.nansum(Omass))
        print('ion fraction =',np.nansum(ion_frac*Omass)/np.nansum(Omass))
    
    med_array.append(np.median(mid_ion_frac))

    print('Now looking at High-mass')
    for s in range(len(lab)-11,len(lab)):
        print('in galaxies',lab[s])

        Omass    = np.loadtxt('/home1/nnsanche/ROMULUS_research/ALLcolumndensities/Omass_z0_17/'+str(lab[s])+'_Omass_3456.np')
        ion_frac = np.loadtxt('/home1/nnsanche/ROMULUS_research/ALLcolumndensities/ion_frac_z0_17/'+str(lab[s])+'_'+ion_labels[i]+'frac_3456.np')

        hi_ion_frac.append(np.nansum(ion_frac*Omass)/np.nansum(Omass))
        print('ion fraction =',np.nansum(ion_frac*Omass)/np.nansum(Omass))

    med_array.append(np.median(hi_ion_frac))
#    quit()
    
    ion_array.append(med_array)

    if i == 0 :
        print('Plotting oi')
        plt.bar(indec,ion_array[i],width=wid,label=ion_labels[i],color=colors[i],edgecolor='Black')
        ion_butt = np.array(ion_array[i]) 
        print('ion butts:',ion_butt)
        plt.text(0.6,0.1,str('%.3f' % ion_array[i][0]),color=colors[i],size=12)
        plt.text(2.6,0.1,str('%.3f' % ion_array[i][1]),color=colors[i],size=12)
        plt.text(4.6,0.1,str('%.3f' % ion_array[i][2]),color=colors[i],size=12)
        t = 0.1

    else :
        print('Plotting other ions')
        plt.bar(indec,ion_array[i],width=wid,label=ion_labels[i],color=colors[i],edgecolor='Black',bottom=ion_butt)
        ion_butt = np.array(ion_butt) + np.array(ion_array[i])
        print('ion butts:',ion_butt)
        plt.text(0.6,0.1+t,str('%.3f' % ion_array[i][0]),color=colors[i],size=12)
        plt.text(2.6,0.1+t,str('%.3f' % ion_array[i][1]),color=colors[i],size=12)
        plt.text(4.6,0.1+t,str('%.3f' % ion_array[i][2]),color=colors[i],size=12)
        t = t + 0.1

plt.text(-0.7,-0.06,r'5 $\times$ 10$^{10}$ - 5 $\times$ 10$^{11}$',size=11)
plt.text(1.3,-0.06,r'5 $\times$ 10$^{11}$ - 2 $\times$ 10$^{12}$ ',size=11)
plt.text(3.3,-0.06,r'2 $\times$ 10$^{12}$ - 2 $\times$ 10$^{13}$ ',size=11)
plt.text(1,-0.12,r'Averaged M$_{halo}$ Range',size=13)
plt.ylabel(r'Oxygen Fraction (CGM <R$_{200}$)',size=13)
plt.ylim(0,1.15)
plt.xticks([])
plt.xlim(-0.75,5.25)
plt.legend(ncol=4)
plt.savefig('OPPEN_fig10_ROMgxys_averaged.pdf')
plt.show()


