""" Creat ensembles
Input:  nc_in     - string.   The nc file you want to copy
        nc_var    - string.   The variable in question, e.g. 'temp'
        number=1  - integer.  The number of ensembles to create
        degree=0.1- float.    The degree of change for the grid box
Output: ens_dict  - dict      Dictionary, keys: ens0, ens1 etc
"""

import numpy as np
import inout as io
from copy import copy


def butterfly(nc_in,nc_var,number=1,degree=0.1):
#x=55, y=25
  data = io.readnetcdf(nc_in,nc_var)[0]
  ens_dict = {}
  ens=np.zeros(data.shape)
  for i in xrange(number):
    n=25+i
    print(n)
    ens=data.copy()
    ens[0,n,55]=data[0,n,55]+degree+0.01*np.random.rand()
    ens_dict[i]=ens.copy()
  return ens_dict

def datafly(data,number=1,degree=0.1):
#x=55, y=25
  ens_dict = {}
  ens=np.zeros(data.shape)
  for i in xrange(number):
    n=25+i
    print(n)
    ens=data.copy()
    ens[0,n,55]=data[0,n,55]+degree+0.01*np.random.rand()
    ens_dict[i]=ens.copy()
  return ens_dict
