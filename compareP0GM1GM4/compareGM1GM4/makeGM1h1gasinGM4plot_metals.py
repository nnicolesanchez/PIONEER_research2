# This script makes some plots                                                                                                    
from scipy.interpolate import griddata
import scipy.interpolate
import matplotlib.pyplot as plt
import numpy as np

new_GM4_gaspos  = np.loadtxt('GM4_GM1matchedgas_posxyz.txt')
new_GM4_gastemp = np.loadtxt('GM4_GM1matchedgas_temp.txt')
new_GM4_gasmetal = np.loadtxt('GM4_GM1matchedgas_metal.txt')

new_GM1_gaspos  = np.loadtxt('GM1matchedgas_posxyz.txt')
new_GM1_gastemp = np.loadtxt('GM1matchedgas_temp.txt')
new_GM1_gasmetal = np.loadtxt('GM1matchedgas_metal.txt')


fig = plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(121) 
ax2 = fig.add_subplot(122) 
ax1.hist(new_GM4_gasmetal,bins=50,range=(-0.1,0.1))
ax2.hist(new_GM1_gasmetal,bins=50,range=(-0.1,0.1))
plt.show()

print('Min gas metallicity',min(new_GM4_gasmetal),'Max gas metallicity',max(new_GM4_gasmetal))
print('# of Gas particles w/ metalfrac < -0.1',len(new_GM4_gasmetal[new_GM4_gasmetal < -0.025]))
print('# of Gas particles w/ metalfrac > 0.1',len(new_GM4_gasmetal[new_GM4_gasmetal > 0.025]))
print('Mean Z',np.mean(new_GM4_gasmetal))

plt.scatter(new_GM4_gaspos[:,0],new_GM4_gaspos[:,1],c=new_GM4_gasmetal,s=10,cmap=plt.cm.get_cmap('jet'),alpha=0.5,vmin=-0.025,vmax=0.025)
plt.title('GM4 Gas Particles (Matched with GM1 Gas)')
plt.xlabel('X [kpc]')
plt.ylabel('Y [kpc]')
plt.ylim(-2000,1200)
plt.xlim(-1900,1800)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Metal Mass Fraction')
plt.savefig('GM4_xy_Z.pdf')
plt.show()

x = new_GM4_gaspos[:,0]
y = new_GM4_gaspos[:,1]
z = np.log10(new_GM4_gasmetal)

xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet,alpha=0.5,vmin=4,vmax=6)
plt.xlabel('X [kpc]')
plt.ylabel('Y [kpc]')
plt.ylim(-2000,1200)
plt.xlim(-1900,1800)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Metal Mass Fration')
plt.title('GM4 Gas Particles (Matched with GM1 Gas)')
plt.savefig('GM4_xy_Zcon.pdf')
plt.show()

#################
##### GM 1 ######
#################

print('Min gas metallicity',min(new_GM1_gasmetal),'Max gas metallicity',max(new_GM1_gasmetal))
print('# of Gas particles w/ Z < -0.25',len(new_GM1_gasmetal[new_GM1_gasmetal < -0.025]))
print('# of Gas particles w/ Z > 0.25',len(new_GM1_gasmetal[new_GM1_gasmetal > 0.025]))
print('Mean Z',np.mean(new_GM1_gasmetal))

plt.scatter(new_GM1_gaspos[:,0],new_GM1_gaspos[:,1],c=new_GM1_gasmetal,s=10,cmap=plt.cm.get_cmap('jet'),alpha=0.5,vmin=-0.025,vmax=0.025)
plt.title('GM1 Gas Particles (Matched with GM4 Gas)')
plt.xlabel('X [kpc]')
plt.ylabel('Y [kpc]')
plt.ylim(-2000,1200)
plt.xlim(-1900,1800)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Metal Mass Fraction')
plt.savefig('GM1_xy_Z.pdf')
plt.show()

x = new_GM1_gaspos[:,0]
y = new_GM1_gaspos[:,1]
z = np.log10(new_GM1_gasmetal)

xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet,alpha=0.5,vmin=4,vmax=6)
plt.xlabel('X [kpc]')
plt.ylabel('Y [kpc]')
plt.ylim(-2000,1200)
plt.xlim(-1900,1800)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Metal Mass Fraction')
plt.title('GM1 Gas Particles (Matched with GM4 Gas)')
plt.savefig('GM1_xy_Zcon.pdf')
plt.show()
