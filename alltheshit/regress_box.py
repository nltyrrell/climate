import numpy as np
import scipy as sp
from scipy import stats
import Scientific.IO.NetCDF as S
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm

"""
@author Nicholas Tyrrell
@brief Regress a timeseries onto a 2D array

Input:   NetCDF file, one 2D variable
Output: 2D array (or NetCDF file?) of 'slope' value (a)

At each gridbox calculate the regression between the surface temperture and a chosen variable, e.g. DLWR

Regress Y onto X, where Y = aX + b.
Y is timeseries, X is 2D array
"""

#txtar = readtext(txt)
#ncdata, ncobj = readnetcdf(ncfile,var)
#reg_map = regmap(txtar,ncdata,statvar,col,lag)
#writenc(outnc,reg_map)

# Regression
"""Maps the regression values (or other values) between a timeseries, e.g. NINO3 and an array, e.g. Hadisst
txt_ts is a text file of the timeseries, and col is the column to use (if there are more than one, otherwise leave blank)
statvar is the statistical variable: slope=0, intercept=1, r_val=2, p_val=3, stderr=4
Lag is in MONTHS 
"""
def regmap(temp,var,stat=0,lag=0):
	t,lat,lon = temp.shape
	reg_map = np.zeros((lat,lon))
	for i in xrange(lat):
		for j in xrange(lon):
			#slope, intercept, r_val, p_val, stderr = stats.linregress(map_ts[:,i,j],txt_ts[:,col])
			#stattup = stats.linregress(map_ts[lag:1706,i,j],txt_ts[0:1706-lag,col])
			stattup = stats.linregress(var[0:t-lag,i,j],temp[lag:t,i,j])
			#stattup = stats.linregress(txt_ts[11+lag:1695,col],map_ts[11:1695-lag,i,j])
			reg_map[i,j] = stattup[stat]
			reg_II = np.nan_to_num(reg_map)
			out_reg = np.ma.masked_outside(reg_II,-1000,1000)
	return out_reg,reg_map,reg_II


## Correlation
#def cormap(txt_ts,map_ts,lag=0):
#  t,lat,lon = map_ts.shape
#  cor_map = np.zeros((lat,lon))
#  for i in xrange(lat):
#    for j in xrange(lon):
#      #slope, intercept, r_val, p_val, stderr = stats.linregress(map_ts[:,i,j],txt_ts[:,col])
#      r_mat = np.corrcoef(map_ts[lag:1706,i,j],txt_ts[0:1706-lag,0])
#      cor_map[i,j] = r_mat[1,0]  
#  out_cor_map = np.nan_to_num(cor_map)
#  return out_cor_map
