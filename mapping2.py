import Scientific.IO.NetCDF as S
import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy as sc
from scipy import stats
from mycmaps import jetwhite
from matplotlib.colors import LogNorm


def linsimplot(var, clf=True):
	plt.ion()
	if clf==True:
		plt.clf()
	plt.plot(var)
	plt.grid()

def simplot(var,clf=True,clim=False):
    """ Creates a pcolormesh plot
    Usage: simplot(array,clim=1.5)
    kw clim gives positive and negative limits"""

    dims = len(var.shape)
    if dims == 4:
            var = var[0,0,::]
    elif dims == 3:
            var = var[0,::]
    else:
            var=var
    plt.ion()
    if clf==True:
            plt.clf()
    if clim!=False:
            plt.pcolormesh(var,cmap=jetwhite(n=256)).set_clim(vmin=-clim,vmax=clim)
    if clim==False:
            plt.pcolormesh(var)
    if clf==True:
            plt.colorbar()
    return

