# This script creates a plot of the CGM mass 
# of different phases:
#        Cool: 10^4   - 10^4.5
#              10^4.5 - 10^5
#        Warm: 10^5   - 10^6
#        Hot:  10^6   - 10^7

#     - Outputs:
#         1. Plot of warm CGM Mass vs radius for 
#            all sims (P0-GM7)
#               - one for cool, c/w, warm, hot


# N. Nicole Sanchez -- September 2017
# Univ. of Wash.    -- Nbody Shop
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np
import sys

plt.rc('font', size=12, family='serif', style='normal', variant='normal', stretch='normal', weight='normal')
plt.rc('xtick', labelsize=12)
plt.rc('xtick.major', size=6, width=1)
plt.rc('lines', lw=2)
plt.rc('axes', lw=1, labelsize=12)

names     = ['P0','GM1','GM4','GM5','GM6','GM7']
colors    = ['DodgerBlue','SteelBlue','FireBrick','IndianRed','Salmon','Orange']
labels    = ['cool','cool/warm','warm','hot']
R_vir      = [262.25,264.5,250,254.7,239.8]
lines     = [':','-.','--','-',':','-.']

#steps1    = ['0454','0972','1739','2554','3195','4096'] # Use 3195 for P0 and 3200 for the rest of the GM
#steps2    = ['0454','0972','1739','2554','3200','4096'] 
#labels2   = ['z = 4','z = 2','z = 1','z = 0.5','z = 0.25','z = 0']

colors2   = ['LightGreen','MediumSeaGreen','Green','DarkSlateGrey']
steps1    = ['0972','1739','2554','4096'] # Use 3195 for P0 and 3200 for the rest of the GM
steps2    = ['0972','1739','2554','4096']
labels2   = ['z = 2','z = 1','z = 0.5','z = 0']


for i in range(len(names)-1):
    fig = plt.figure(figsize=(8, 5))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
#    ax3 = fig.add_subplot(313)

    if (i == 0):
        steps = steps1
    else:
        steps = steps2

    for j in range(len(steps)):
        radii = np.loadtxt(names[i]+'_radius.txt')
        total = np.loadtxt(names[i]+'_totCGMmass_'+steps[j]+'.txt')
        warm  = np.loadtxt(names[i]+'_warmCGMmass_'+steps[j]+'.txt')
        warm1 = np.loadtxt(names[i]+'_warm_half1CGMmass_'+steps[j]+'.txt')
        warm2 = np.loadtxt(names[i]+'_warm_half2CGMmass_'+steps[j]+'.txt')

        ax1.plot(radii[:-1],np.log10(warm1),color=colors2[j],linestyle=lines[j],label=str(labels2[j]))
        ax2.plot(radii[:-1],np.log10(warm2),color=colors2[j],linestyle=lines[j],label=str(labels2[j]))
        ax1.plot([R_vir[i],R_vir[i]],[4,10],color='grey',linestyle='--')
        ax2.plot([R_vir[i],R_vir[i]],[4,10],color='grey',linestyle='--')
#        ax3.plot(radii[:-1],np.log10(warm),color=colors2[j],linestyle=lines[j],label=str(labels2[j]))


#    ax3.set_ylabel(r'T = 10$^{5}$ - 10$^{6}$ K')
    ax1.set_ylabel(r'T = 10$^{5-5.7}$ K')
    ax2.set_ylabel(r'T = 10$^{5.7-6}$ K')
    ax1.set_title(str(names[i])+' Warm CGM Mass Profile [M$_{\odot}$]')
    ax1.set_ylim(5,9)
    ax2.set_ylim(5,9)
    ax1.set_xlim(-10,270)
    ax2.set_xlim(-10,270)
    plt.xlabel('Radius [kpc]')
    plt.subplots_adjust(hspace = .001)
    ax1.legend(loc='right')#bbox_to_anchor=(.95, 1.0))
    plt.savefig(names[i]+'_warmCGMmassvsradius_split.pdf')
    plt.show()
    plt.close()
#    quit()
quit()

fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

for i in range(len(names)):
    for j in range(len(steps)):
        radii = np.loadtxt(names[i]+'_radius.txt')
        total = np.loadtxt(names[i]+'_totCGMmass_'+steps[j]+'.txt')
        warm  = np.loadtxt(names[i]+'_warmCGMmass_'+steps[j]+'.txt')
        warm1 = np.loadtxt(names[i]+'_warm_half1CGMmass_'+steps[j]+'.txt')
        warm2 = np.loadtxt(names[i]+'_warm_half2CGMmass_'+steps[j]+'.txt')
        

        ax1.plot(radii[:-1],np.log10(warm),color=colors[i],linestyle=lines[i],label=names[i])
        ax2.plot(radii[:-1],np.log10(warm1),color=colors[i],linestyle=lines[i],label=names[i])
        ax3.plot(radii[:-1],np.log10(warm2),color=colors[i],linestyle=lines[i],label=names[i])

plt.ylabel(r'CGM Cumulative Mass [M$_{\odot}$]')
plt.xlabel('Radius [kpc]')
plt.legend()
plt.savefig('ALLwarmCGMmassvsradius.pdf')
plt.show()
    
