# This script reads in the data created by find_GM1h1gasinGM4.py
# which prints out the matched gas particles between GM1 & GM4
# and creates:
#         - X vs Y plots for GM4 and GM1
#         ** Colors corresponding to tempertaure
#         ** Second plots is a attempted contour plot
#               Fix later because still looks messy


# N. Nicole Sanchez -- June 23 2017
# Univ. of Wash.    -- Nbody Shop
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

################
##### GM4 ######
################

print('Min gas temp',min(new_GM4_gastemp),'Max gas temp',max(new_GM4_gastemp))
print('# of Gas particles w/ T < 10^4 K',len(new_GM4_gastemp[new_GM4_gastemp < 10**4]))
print('# of Gas particles w/ T > 10^6 K',len(new_GM4_gastemp[new_GM4_gastemp > 10**6]))
print('% of Gas particles left out of range',str.format('{0:0.3f}',float(len(new_GM4_gastemp[new_GM4_gastemp < 10**4])+len(new_GM4_gastemp[new_GM4_gastemp > 10**6]))/len(new_GM4_gastemp)*100),'%')

#fig = plt.figure(figsize=(15, 5))
#ax1 = fig.add_subplot(121)
#ax2 = fig.add_subplot(122)

plt.scatter(new_GM4_gaspos[:,0],new_GM4_gaspos[:,1],c=np.log10(new_GM4_gastemp),s=10,cmap=plt.cm.get_cmap('jet'),alpha=0.5,vmin=4,vmax=6)
plt.xlabel('X [kpc]')
plt.ylabel('Y [kpc]')
plt.ylim(-2000,1200)
plt.xlim(-1900,1800)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Temperature [log(K)]')
plt.savefig('GM4_xy_T.pdf')
plt.show()

x = new_GM4_gaspos[:,0]
y = new_GM4_gaspos[:,1]
z = np.log10(new_GM4_gastemp)

xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

# Contour the gridded data, plotting dots at the randomly spaced data points.
#CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet,alpha=0.5,vmin=4,vmax=6)
plt.xlabel('X [kpc]')
plt.ylabel('Y [kpc]')
plt.ylim(-2000,1200)
plt.xlim(-1900,1800)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Temperature [log(K)]')
plt.title('GM4 Gas Particles (Matched with GM1 Gas)')
plt.savefig('GM4_xy_Tcon.pdf')
plt.show()

################
##### GM1 ######
################

print('Min gas temp',min(new_GM1_gastemp),'Max gas temp',max(new_GM1_gastemp))
print('# of Gas particles w/ T < 10^4 K',len(new_GM1_gastemp[new_GM1_gastemp < 10**4]))
print('# of Gas particles w/ T > 10^6 K',len(new_GM1_gastemp[new_GM1_gastemp > 10**6]))
print('% of Gas particles left out of range',str.format('{0:0.3f}',float(len(new_GM1_gastemp[new_GM1_gastemp < 10**4])+len(new_GM1_gastemp[new_GM1_gastemp > 10**6]))/len(new_GM1_gastemp)*100),'%')

plt.scatter(new_GM1_gaspos[:,0],new_GM1_gaspos[:,1],c=np.log10(new_GM1_gastemp),s=10,cmap=plt.cm.get_cmap('jet'),alpha=0.5,vmin=4,vmax=6)
plt.title('GM1 Gas Particles (Matched with GM4 Gas)')
plt.xlabel('X [kpc]')
plt.ylabel('Y [kpc]')
plt.ylim(-2000,1200)
plt.xlim(-1900,1800)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Temperature [log(K)]')
plt.savefig('GM1_xy_T.pdf')
plt.show()


x = new_GM1_gaspos[:,0]
y = new_GM1_gaspos[:,1]
z = np.log10(new_GM1_gastemp)

xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet,alpha=0.5,vmin=4,vmax=6)
plt.xlabel('X [kpc]')
plt.ylabel('Y [kpc]')
plt.ylim(-2000,1200)
plt.xlim(-1900,1800)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Temperature [log(K)]')
plt.title('GM1 Gas Particles (Matched with GM4 Gas)')
plt.savefig('GM1_xy_Tcon.pdf')
plt.show()
