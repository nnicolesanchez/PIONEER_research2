# This script plots the star formation history for the GM                                           
# suite for h243 given OUTPUTS from sfhistory.py  
import matplotlib.pyplot as plt
import numpy as np
import pynbody

plt.rc('font', size=12, family='serif', style='normal', variant='normal', stretch='normal', weight='normal')
plt.rc('xtick', labelsize=12)
plt.rc('xtick.major', size=6, width=1)
plt.rc('lines', lw=2)
plt.rc('axes', lw=1, labelsize=12)

data =['P0_sfhistory_bins.txt','GM1_sfhistory_bins.txt','GM4_sfhistory_bins.txt','GM5_sfhistory_bins.txt','GM6_sfhistory_bins.txt','GM7_sfhistory_bins.txt']
labels = ['P0','GM1','GM4','GM5','GM6','GM7']
colors = ['DodgerBlue','SteelBlue','FireBrick','IndianRed','Salmon','Orange']

ticks_P0  = [3,9,12,15,21]
ticks_GM1 = [5,18,31,44,69]
ticks_GM456  = [5,14,21,28,42]
ticks     = [ticks_P0,ticks_GM1,ticks_GM456,ticks_GM456,ticks_GM456,ticks_GM456]

fig = plt.figure()
ax1 = fig.add_subplot(111)
#ax2 = ax1.twiny()

for i in range(len(data)):
    dt = np.transpose(np.loadtxt(data[i]))
    z = np.array(np.loadtxt('../compareALL/'+labels[i]+'_redshifts.txt'))
    SFR = dt[0]
    t = dt[1]

    plt.plot(t,SFR,color=colors[i],label=labels[i])#,drawstyle='steps')

    print(t)
    print(ticks[i][2])
    print(t[ticks[i][2]])
    print(z[ticks[i][2]])
    new_tick_locations = [t[ticks[i][0]],t[ticks[i][1]],t[ticks[i][2]],t[ticks[i][3]],t[ticks[i][4]]]
    new_tick_labels = ["%.0f" % z[ticks[i][0]],"%.0f" % z[ticks[i][1]],"%.1f" % z[ticks[i][2]],"%.1f" % z[ticks[i][3]],"%.1f" % np.abs(z[ticks[i][4]])]

#ax2.set_xticks(new_tick_locations)
#ax2.set_xticklabels(new_tick_labels)
#ax2.set_xlabel('z')

ax1.set_ylim(0,6)
ax1.set_xlim(5,6.5)
#ax2.set_xlim(5,6.5)
ax1.set_ylabel(r'SFR [M$_{\odot}$ yr$^{-1}$]')
ax1.set_xlabel('Time [Gyr]')
plt.legend()
plt.savefig('allSFH_zoomin.pdf')
plt.show()
