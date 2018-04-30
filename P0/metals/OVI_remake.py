# In-the-works test remake of OVI_frompontzen.py
# Switching to a new snapshot because 000832 is being weird?
import matplotlib.pyplot as plt
import warnings
import pynbody
import numpy as np


#plt.ion()
# This @pynbody systems basically creates the array you need assuming it
# is dependent on some other array you already have
# a. allows code to assume this secondary array exists
# b. keeps these array updated (since they are always based on initial array)
@pynbody.derived_array
def rhoOVI(f):
# This returns the array: OVI density
# looks like gas density * OVI ion fractions * oxygen mass fraction
    ovi = pynbody.analysis.ionfrac.calculate(f,ion='ovi',mode='new')
    # so ovi is an array of the calculate ion fractions in this snapshot
    # designated f and loaded AFTER since this array will exist is f does
    return f.gas['rho']*ovi*f.gas['OxMassFrac']

# Load in the snapshot
#f = pynbody.load("/u/apontzen/scratch/pioneer50h128.1536gst1.bwK1/pioneer50h128.1536gst1.bwK1.000384") 
f = pynbody.load("/u/apontzen/scratch/pioneer50h128.1536gst1.bwK1/pioneer50h128.1536gst1.bwK1.000128")
#f = pynbody.load("/u/apontzen/scratch/pioneer50h128.1536gst1.bwK1/pioneer50h128.1536gst1.bwK1.001024")

# This appears to center the halo around the stars?
pynbody.analysis.halo.center(f.star)
print(f.properties)
print(len(f.star))

# Translates everything into physical units
f.physical_units()

# Plots stuff (units are not meaningful to me)
pynbody.plot.sph.image(f.gas,qty='rhoOVI',width='500 kpc',units="16 m_p cm^-2",vmin=1e13,vmax=1e15)#,cmap='magma')
plt.savefig('rho_ovi_image.png')
plt.show()
