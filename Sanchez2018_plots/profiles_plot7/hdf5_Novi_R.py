from pynbody.analysis import profile
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import seaborn as sns
import numpy as np
import pynbody
import h5py
import sys

#from ..array import SimArray
from pynbody import config
import logging
logger = logging.getLogger('pynbody.analysis.ionfrac')

#from scipy.interpolate import interp3d
from scipy.interpolate import RegularGridInterpolator as rgi

from pynbody.analysis.interpolate import _interpolate3d

def interpolate3d(x, y, z, x_vals, y_vals, z_vals, vals):
    """
    Interpolate on a 3D regular grid. 
    Yields results identical to scipy.interpolate.interpn. 

    Input
    -----

    x,y,z : points where the interpolation will be performed

    x_vals, y_vals, z_vals : xyz values of the reference grid

    vals : grid values
    """

    # cast x_vals, y_vals and z_vals to float64

    x_vals = x_vals.astype(np.float64)
    y_vals = y_vals.astype(np.float64)
    z_vals = z_vals.astype(np.float64)
    vals = vals.astype(np.float64)

    result_array = np.empty(len(x), dtype=np.float64)

    _interpolate3d.interpolate3d(len(x),
                                 x, y, z,
                                 len(x_vals), x_vals,
                                 len(y_vals), y_vals,
                                 len(z_vals), z_vals,
                                 vals,
                                 result_array)

    return result_array

def N_OVI(f):
    ovi = pynbody.analysis.ionfrac_edit.calculate(f,ion='ovi',mode='new')
    m_p = 1.6726 * 10**-24 #g
    #print(ovi)
    #print(f.gas['OxMassFrac'])
    #print(f.gas['rho'].in_units('g cm**-3')*ovi*f.gas['OxMassFrac']/(16*m_p))
    return f.gas['rho'].in_units('g cm**-3')*ovi*f.gas['OxMassFrac']/(16*m_p)


def hdf5_ion_frac(sim, ion):
    if ion == 'oi':
        print('Loading OI')
        nion = 1
    elif ion == 'oii':
        print('Loading OII')
        nion = 2
    elif ion == 'oiii':
        print('Loading OIII')
        nion = 3
    elif ion == 'oiv':
        print('Loading OIV')
        nion = 4
    elif ion == 'ov':
        print('Loading OV')
        nion = 5
    elif ion == 'ovi':
        print('Loading OVI')
        nion = 6
    elif ion == 'ovii':
        print('Loading OVII')
        nion = 7
    elif ion == 'oviii':
        print('Loading OVIII')
        nion = 8
    else:
        print('Specified ion incompatible; Try again.')
        

    iffile = '/home1/nnsanche/hm2012_hr.h5'
    ifs = h5py.File(iffile)
    
    x_vals = ifs['O'].attrs['Parameter2']   # Redshifts
    y_vals = ifs['O'].attrs['Temperature']  # Temperatures
    z_vals = ifs['O'].attrs['Parameter1']   # Densities

    vals = ifs['O'][nion]
    x = np.zeros(len(sim.gas))
    x[:] = sim.properties['z']
    x = x
    y = np.log10(sim.gas['temp']).view(np.ndarray)
    z = np.log10(sim.gas['rho'].in_units('m_p cm^-3')).view(np.ndarray)
 
    n = len(sim.gas)
    n_z_vals = len(z_vals)
    n_y_vals = len(y_vals)
    n_x_vals = len(x_vals)
    
    # get values off grid to minmax                 
    x[np.where(x < np.min(x_vals))] = np.min(x_vals)
    x[np.where(x > np.max(x_vals))] = np.max(x_vals)
    y[np.where(y < np.min(y_vals))] = np.min(y_vals)
    y[np.where(y > np.max(y_vals))] = np.max(y_vals)
    z[np.where(z < np.min(z_vals))] = np.min(z_vals)
    z[np.where(z > np.max(z_vals))] = np.max(z_vals)
    
    # interpolate                              
#    logger.info("Interpolation %s values" % ion)
    result_array = interpolate3d(x, y, z, x_vals, y_vals, z_vals, vals)
#    my_interpolating_function = rgi((x,y,z), vals)
#    result_array = my_interpolating_function(array([x_vals,y_vals,z_vals]).T)

    return 10 ** result_array





