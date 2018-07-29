import matplotlib
matplotlib.use('Agg')


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
# See link for cookbook and details
# http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay
# Thank you for introducing the Savitzky-Golay filter! So basically this is just like a regular "Moving average" filter, but instead of just calculating the average, a polynomial (usually 2nd or 4th order) fit is made for every point, and only the "middle" point is chosen. Since 2nd (or 4th) order information is concerned at every point, the bias introduced in "moving average" approach at local maxima or minima, is circumvented. Really elegant
    import numpy as np
    from math import factorial
  
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def smooth(y, box_pts):
    # moving average smoothing
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


import tangos
tangos.all_simulations()
import pylab as p
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

P0_halo = tangos.get_halo("pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.004096/halo_1")
P0_SFR = P0_halo["SFR_histogram"]
P0_BH_accrate = P0_halo.calculate('BH.BH_mdot_histogram')
#P0_BH_accrate_hat = savitzky_golay(P0_BH_accrate, 51, 3)
P0_BH_accrate_hat = smooth(P0_BH_accrate, 30)
P0_SFR_property_object = P0_halo.get_objects("SFR_histogram")[0]
P0_SFR_time_bins = P0_SFR_property_object.x_values()

GM1_halo = tangos.get_halo("pioneer50h243GM1.1536gst1bwK1BH_no3072/pioneer50h243GM1.1536gst1bwK1BH.004096/halo_1")
GM1_SFR = GM1_halo["SFR_histogram"]
GM1_BH_accrate = GM1_halo.calculate('BH.BH_mdot_histogram')
#GM1_BH_accrate_hat = savitzky_golay(GM1_BH_accrate, 51, 3)
GM1_BH_accrate_hat = smooth(GM1_BH_accrate, 30)
GM1_SFR_property_object = GM1_halo.get_objects("SFR_histogram")[0]
GM1_SFR_time_bins = GM1_SFR_property_object.x_values()

GM7_halo = tangos.get_halo("pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.003968/halo_1")
GM7_SFR = GM7_halo["SFR_histogram"]
GM7_BH_accrate = GM7_halo.calculate('BH.BH_mdot_histogram')
#GM7_BH_accrate_hat = savitzky_golay(GM7_BH_accrate, 51, 3)
GM7_BH_accrate_hat = smooth(GM7_BH_accrate, 30)
GM7_SFR_property_object = GM7_halo.get_objects("SFR_histogram")[0]
GM7_SFR_time_bins = GM7_SFR_property_object.x_values()

GM4_halo = tangos.get_halo("pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.004096/halo_1")
GM4_SFR = GM4_halo["SFR_histogram"]
GM4_BH_accrate = GM4_halo.calculate('BH.BH_mdot_histogram')
#GM4_BH_accrate_hat = savitzky_golay(GM4_BH_accrate, 51, 3)
GM4_BH_accrate_hat = smooth(GM4_BH_accrate, 30)
GM4_SFR_property_object = GM4_halo.get_objects("SFR_histogram")[0]
GM4_SFR_time_bins = GM4_SFR_property_object.x_values()


plt.plot(P0_SFR_time_bins, np.log10(P0_BH_accrate_hat),color='DodgerBlue',label='P0')
plt.plot(GM1_SFR_time_bins, np.log10(GM1_BH_accrate_hat),color='SteelBlue',label='GM1')
plt.plot(GM7_SFR_time_bins, np.log10(GM7_BH_accrate_hat),color='FireBrick',label='GM2')
plt.plot(GM4_SFR_time_bins, np.log10(GM4_BH_accrate_hat),color='Salmon',label='GM3')
plt.xlabel("Age/Gyr")
plt.ylabel("BH accretion rate/$M_{\odot}\,yr^{-1}$")
plt.ylim(-6,-1)
plt.legend()
plt.savefig('ALL_bhaccr_age.pdf')
plt.show()
plt.clf()

P0_BH = P0_halo['BH'][0]
GM1_BH = GM1_halo['BH'][0]
GM7_BH = GM7_halo['BH'][0]
GM4_BH = GM4_halo['BH'][0]

print(GM7_BH.keys())
plt.plot(P0_BH.calculate_for_progenitors('t()')[0][:-1],np.log10(P0_BH.calculate_for_progenitors('BH_mass')[0]),color='DodgerBlue',label='P0')
plt.plot(GM1_BH.calculate_for_progenitors('t()')[0][:-3],np.log10(GM1_BH.calculate_for_progenitors('BH_mass')[0]),color='SteelBlue',label='GM1')
plt.plot(GM7_BH.calculate_for_progenitors('t()')[0],np.log10(GM7_BH.calculate_for_progenitors('BH_mass')[0]),color='FireBrick',label='GM2')
plt.plot(GM4_BH.calculate_for_progenitors('t()')[0][:-2],np.log10(GM4_BH.calculate_for_progenitors('BH_mass')[0]),color='Salmon',label='GM3')
plt.ylabel(r"log(M$_{BH}$/M$_{\odot}$)")
plt.xlabel("Age/Gyr")
plt.legend()
plt.savefig('ALL_bhmass_age.pdf')
plt.show()
