# This script creates a plot of the CGM mass 
# of different phases:
#        Cool: 10^4   - 10^4.5
#              10^4.5 - 10^5
#        Warm: 10^5   - 10^5.7
#              10^5.7 - 10^6
#        Hot:  10^6   - 10^7

#     - Outputs:
#         1. Plot of CGM Mass vs redshift/time for 
#            all sims (P0-GM7)
#               - one for cool, c/w, warm, hot
#         2. Plot of CGM Mass vs redshift for EACH sim
#               - one of cool/c/w/warm/hot for P0, then GM1 etc
#         3. Mean CGM metal fraction vs redshift/time
#               - for each sim (divided by phase, like #2)

# N. Nicole Sanchez -- August 2017
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
lines     = ['-','--','-.',':']
colors2   = ['lightcoral','indianred','FireBrick','DarkRed']
ticks_P0  = [3,9,12,15,21]
ticks_GM1 = [5,18,31,44,69]
ticks_GM456  = [5,14,21,28,42]
ticks     = [ticks_P0,ticks_GM1,ticks_GM456,ticks_GM456,ticks_GM456]

# Load in the data and create phase arrays
cool  = []
cw    = []
warm1 = []
warm2 = []
hot   = []
t     = []
z     = []
stars = []

for j in range(4):#len(names)):
    cool.append(np.loadtxt(names[j]+'_phaseCGM_mass_cool.txt'))
    cw.append(np.loadtxt(names[j]+'_phaseCGM_mass_cool_warm.txt'))
    warm1.append(np.loadtxt(names[j]+'_phaseCGM_mass_warm_half1.txt'))
    warm2.append(np.loadtxt(names[j]+'_phaseCGM_mass_warm_half2.txt'))
    hot.append(np.loadtxt(names[j]+'_phaseCGM_mass_hot.txt'))
    t.append(np.loadtxt(names[j]+'_time_Gyr.txt'))
    z.append(np.loadtxt(names[j]+'_redshifts.txt'))
    stars.append(np.loadtxt(names[j]+'_totstellarmass.txt'))
cool45 = cool + cw



# CGM Mass vs Redshift/Time per Simulation
for i in range(4):#len(names)):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()

    print(i)
    print(t)
    ax1.plot(t[i],np.log10(hot[i]),color=colors2[3],label=r'CGM 10$^6$ - 10$^7$ K',linestyle=lines[3])
    ax1.plot(t[i],np.log10(warm2[i]),color=colors2[2],label=r'CGM 10$^{5.7}$ - 10$^6$ K',linestyle=lines[2])
    ax1.plot(t[i],np.log10(warm1[i]),color=colors2[1],label=r'CGM 10$^5$ - 10$^{5.7}$ K',linestyle=lines[1])
    ax1.plot(t[i],np.log10(cool45[i]),color=colors2[0],label=r'CGM 10$^4$ - 10$^5$ K',linestyle=lines[0])
    ax1.plot(t[i],np.log10(stars[i]),color=colors[1],label='Stellar Mass',linestyle=lines[1])
    ax1.set_ylabel('log(CGM Mass by Phase)')
    ax1.set_xlabel('Time [Gyr]')

    new_tick_locations = [t[i][ticks[i][0]],t[i][ticks[i][1]],t[i][ticks[i][2]],t[i][ticks[i][3]],t[i][ticks[i][4]]]
    new_tick_labels = ["%.0f" % z[i][ticks[i][0]],"%.0f" % z[i][ticks[i][1]],"%.0f" % z[i][ticks[i][2]],"%.1f" % z[i][ticks[i][3]],"%.0f" % np.abs(z[i][ticks[i][4]])]
#    new_tick_labels = ["%.3f" % z[i][ticks[i][0]],"%.3f" % z[i][ticks[i][1]],"%.3f" % z[i][ticks[i][2]],"%.3f" % z[i][ticks[i][3]],"%.3f" % np.abs(z[i][ticks[i][4]])]
    print(new_tick_locations)
    print(new_tick_labels)

    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(new_tick_labels)
    ax2.set_xlabel('z')

# FULL PLOT
    ax1.set_xlim(-1, 15)
    ax2.set_xlim(-1, 15)
    ax1.set_ylim(5.25,11.25)
    plt.text(-0.5,10.75,names[i],size=12)
    plt.savefig(str(names[i])+'_CGMmassbyphase_redshift.pdf')

# ZOOM IN
#    ax1.set_xlim(3.25, 9) 
#    ax2.set_xlim(3.25, 9) 
#    ax1.set_ylim(5.25,11.25)
#    plt.text(3.5,10.75,names[i],size=12)
#    plt.savefig(str(names[i])+'zoomin_CGMmassbyphase_redshift.pdf')

    ax1.legend(loc=4)
    plt.show()
    plt.close()

quit()
# COOL CGM Mass (for all gxys) vs Redshift
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

ax1.plot(t[0],np.log10(cool45[0]),color=colors[0],label=names[0],linestyle=lines[3])
ax1.plot(t[1],np.log10(cool45[1]),color=colors[1],label=names[1],linestyle=lines[2])
ax1.plot(t[2],np.log10(cool45[2]),color=colors[2],label=names[2],linestyle=lines[1])
ax1.plot(t[3],np.log10(cool45[3]),color=colors[3],label=names[3],linestyle=lines[0])
ax1.plot(t[4],np.log10(cool45[4]),color=colors[4],label=names[4],linestyle=lines[3])
#ax1.plot(t[5],np.log10(cool[5]),color=colors[5],label=names[5],linestyle=lines[2])
ax1.set_ylabel(r'log(Cool CGM Mass) [T$^4$ - T$^5$]')
ax1.set_xlabel('Time [Gyr]')
ax1.legend()
plt.savefig('ALL_coolCGMmass_redshift.pdf')
plt.show()
plt.close()


# WARM CGM Mass (for all gxys) vs Redshift
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

ax1.plot(t[0],np.log10(warm1[0]),color=colors[0],label=names[0],linestyle=lines[3])
ax1.plot(t[1],np.log10(warm1[1]),color=colors[1],label=names[1],linestyle=lines[2])
ax1.plot(t[2],np.log10(warm1[2]),color=colors[2],label=names[2],linestyle=lines[1])
ax1.plot(t[3],np.log10(warm1[3]),color=colors[3],label=names[3],linestyle=lines[0])
ax1.plot(t[4],np.log10(warm1[4]),color=colors[4],label=names[4],linestyle=lines[3])
#ax1.plot(t[5],np.log10(warm[5]),color=colors[5],label=names[5],linestyle=lines[2])
ax1.set_ylabel(r'log(Warm CGM Mass) [T$^5$ - T$^{5.7}$]')
ax1.set_xlabel('Time [Gyr]')
ax1.legend()
plt.savefig('ALL_warm1CGMmass_redshift.pdf')
plt.show()

# WARM CGM Mass (for all gxys) vs Redshift
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

ax1.plot(t[0],np.log10(warm2[0]),color=colors[0],label=names[0],linestyle=lines[3])
ax1.plot(t[1],np.log10(warm2[1]),color=colors[1],label=names[1],linestyle=lines[2])
ax1.plot(t[2],np.log10(warm2[2]),color=colors[2],label=names[2],linestyle=lines[1])
ax1.plot(t[3],np.log10(warm2[3]),color=colors[3],label=names[3],linestyle=lines[0])
ax1.plot(t[4],np.log10(warm2[4]),color=colors[4],label=names[4],linestyle=lines[3])
#ax1.plot(t[5],np.log10(warm[5]),color=colors[5],label=names[5],linestyle=lines[2])
ax1.set_ylabel(r'log(Warm CGM Mass) [T$^{5.7}$ - T$^6$]')
ax1.set_xlabel('Time [Gyr]')
ax1.legend()
plt.savefig('ALL_warm2CGMmass_redshift.pdf')
plt.show()

# HOT CGM Mass (for all gxys) vs Redshift
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

ax1.plot(t[0],np.log10(hot[0]),color=colors[0],label=names[0],linestyle=lines[3])
ax1.plot(t[1],np.log10(hot[1]),color=colors[1],label=names[1],linestyle=lines[2])
ax1.plot(t[2],np.log10(hot[2]),color=colors[2],label=names[2],linestyle=lines[1])
ax1.plot(t[3],np.log10(hot[3]),color=colors[3],label=names[3],linestyle=lines[0])
ax1.plot(t[4],np.log10(hot[4]),color=colors[4],label=names[4],linestyle=lines[3])
#ax1.plot(t[5],np.log10(hot[5]),color=colors[5],label=names[5],linestyle=lines[2])
ax1.set_ylabel(r'log(Hot CGM Mass) [T$^5$ - T$^6$]')
ax1.set_xlabel('Time [Gyr]')
ax1.legend()
plt.savefig('ALL_hotCGMmass_redshift.pdf')
plt.show()


