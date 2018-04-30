# This script creates colordensity (phase diagrams) for 
# ChaNGa Nbody Simulations


# N. Nicole Sanchez -- July 6, 2017
# U. W. Seattle





X,Y=np.meshgrid(xrange,yrange)
H = np.log10(H)
masked_array = np.ma.array(H, mask=np.isnan(H))  # mask out all nan, i.e. log10(0.0)
cax = (ax2dhist.pcolormesh(X,Y,masked_array, cmap=cmap, norm=LogNorm(vmin=1,vmax=8)))
