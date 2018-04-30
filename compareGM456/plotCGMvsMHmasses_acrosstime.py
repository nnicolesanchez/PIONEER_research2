import matplotlib.pyplot as plt
import pynbody
import numpy as np



GM4_MH = np.loadtxt('GM4_MH_mass.txt')
GM4_coolCGM = np.loadtxt('GM4_coolCGM_mass.txt')
GM5_MH = np.loadtxt('GM5_MH_mass.txt')
GM5_coolCGM = np.loadtxt('GM5_coolCGM_mass.txt')
GM6_MH = np.loadtxt('GM6_MH_mass.txt')
GM6_coolCGM = np.loadtxt('GM6_coolCGM_mass.txt')

names  = ['GM4','GM5','GM6']
colors = ['FireBrick','IndianRed','Salmon']

#print(GM4_MH,GM5_MH,GM6_MH)

plt.plot(GM4_coolCGM,GM4_MH,label=names[0],marker='*',color=colors[0])
plt.plot(GM5_coolCGM,GM5_MH,label=names[1],marker='*',color=colors[1])
plt.plot(GM6_coolCGM,GM6_MH,label=names[2],marker='*',color=colors[2])

plt.text(0.05*10**10,0.5*10**11,'z ~ 10')
plt.text(0.05*10**10,8*10**11,'z = 0')

plt.xlabel(r'Total CGM Mass [M$_{sol}$] between 10$^4$-10$^5$ K')
plt.ylabel('Main Halo Mass [M$_{sol}$]')
plt.ylim(0,8.5*10**11)
plt.xlim(0,1.3*10**10)
plt.legend()
plt.savefig('GM456_MHmassvsCGMcoolgasmass_acrosstime.pdf')
plt.show()
