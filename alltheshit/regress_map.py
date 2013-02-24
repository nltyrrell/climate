import numpy as np
import scipy as sp
from scipy import stats
import Scientific.IO.NetCDF as S
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm

"""
@author Nicholas Tyrrell
@brief Regress a timeseries onto a 2D array

Input:  Timeseries.txt, e.g. NINO3
        NetCDF file, one 2D variable
Output: 2D array (or NetCDF file?) of 'slope' value (a)

Regress Y onto X, where Y = aX + b.
Y is timeseries, X is 2D array
"""

def algorhythm_reg(txt,ncfile,var,vmin=-1,vmax=1,statvar=0,col=0,lag=0,title='title',outfile='test.eps', outnc='test.nc'):
  """ Example: textar = algorythm('nino3_had_rsc.txt','ncfiles/had.rsc.nc','sst',-1,1,0,0)
  """
  txtar = readtext(txt)
  ncdata, ncobj = readnetcdf(ncfile,var)
  reg_map = regmap(txtar,ncdata,statvar,col,lag)
  #cor_map = cormap(txtar,ncdata,lag)
  #makemap(ncobj,reg_map,plt.cm.RdBu_r,title,outfile,vmin,vmax)
  #makemap(ncobj,cor_map,plt.cm.RdBu_r,title,'cormap.eps',vmin,vmax)
  writenc(outnc,reg_map)
  return reg_map, ncdata,ncobj

# Import timeseries
def readtext(txt_file):
  # Function reads in correctly formatted text file of SST anomalies, outputs numpy array
  fileobj = open(txt_file, 'r')
  outputstr = fileobj.readlines()
  fileobj.close()
  textarray = np.array([[float(n) for n in line.split()] for line in outputstr])
  return textarray


# Import map
def readnetcdf(nc_file,variable):
  """ Function the reads in an nc file and outputs 
  a numpy array (I think)
  """
  nc_fileobj = S.NetCDFFile(nc_file,mode='r')
  nc_data = nc_fileobj.variables[variable].getValue()
  return nc_data, nc_fileobj

# Correlation
def cormap(txt_ts,map_ts,lag=0):
  t,lat,lon = map_ts.shape
  cor_map = np.zeros((lat,lon))
  for i in xrange(lat):
    for j in xrange(lon):
      #slope, intercept, r_val, p_val, stderr = stats.linregress(map_ts[:,i,j],txt_ts[:,col])
      r_mat = np.corrcoef(map_ts[lag:1706,i,j],txt_ts[0:1706-lag,0])
      cor_map[i,j] = r_mat[1,0]  
  out_cor_map = np.nan_to_num(cor_map)
  return out_cor_map

# Regression
"""Maps the regression values (or other values) between a timeseries, e.g. NINO3 and an array, e.g. Hadisst
txt_ts is a text file of the timeseries, and col is the column to use (if there are more than one, otherwise leave blank)
statvar is the statistical variable: slope=0, intercept=1, r_val=2, p_val=3, stderr=4
Lag is in MONTHS 
"""
def regmap(txt_ts,map_ts,statvar=0,col=0,lag=0):
  t,lat,lon = map_ts.shape
  reg_map = np.zeros((lat,lon))
  for i in xrange(lat):
    for j in xrange(lon):
      #slope, intercept, r_val, p_val, stderr = stats.linregress(map_ts[:,i,j],txt_ts[:,col])
      #stattup = stats.linregress(map_ts[lag:1706,i,j],txt_ts[0:1706-lag,col])
      stattup = stats.linregress(txt_ts[0:1706-lag,col],map_ts[lag:1706,i,j])
      #stattup = stats.linregress(txt_ts[11+lag:1695,col],map_ts[11:1695-lag,i,j])
      reg_map[i,j] = stattup[statvar]
  reg_II = np.nan_to_num(reg_map)
  out_reg = np.ma.masked_outside(reg_II,-1000,1000)
  return out_reg,reg_map,reg_II


# Export result map
def writenc(nc_fileout,ancil):
  """ Use regmap, or similar, and add anomaly to a 12 month climatology """
  nc_out = S.NetCDFFile(nc_fileout,mode='a')
  var=nc_out.variables['temp']
  data_out=var.getValue()
  anc_exp = np.repeat(np.expand_dims(ancil,axis=0),12,0)
  ancil2=np.expand_dims(anc_exp,axis=1)
  data_out[:,:,:,:]=ancil2[:,:,:,:]
  var[:]=data_out
  nc_out.close()
  return ancil2


# Display result map

def makemap(nc_obj,mapvar,color_code='plt.cm.RdBu_r',title='title',outfile='test.eps',v_min=None,v_max=None, Show=False):
  """
  Map a climate variable using basemap and matplotlib
  for color_code, need to do plt.cm.hot_r etc
  """
  
  if v_min == None:
    v_min=mapvar.mean()-3*mapvar.std()
    print v_min
  elif v_max == None:
    v_max=mapvar.mean()+3*mapvar.std()
    print v_max

  lon=nc_obj.variables['lon'].getValue()
  lon_units = nc_obj.variables['lon'].units
  lat = nc_obj.variables['lat'].getValue()
  lat_units = nc_obj.variables['lat'].units
  #- Make 2-D longitude and latitude arrays:
  [lon2d, lat2d] = np.meshgrid(lon, lat)

  #- Set up map:
  plt.clf()
  mapproj = bm.Basemap(projection='cyl', llcrnrlat=-90.0, llcrnrlon=0.0, urcrnrlat=90.0, urcrnrlon=356.5)
  mapproj.drawcoastlines()
  mapproj.drawparallels(np.array([-90, -45, 0, 45, 90]), labels=[1,0,0,0])
  mapproj.drawmeridians(np.array([0, 90, 180, 270, 360]), labels=[0,0,0,1])
  lonall, latall = mapproj(lon2d, lat2d)
  plt.ion()
  themap = plt.pcolormesh(lonall, latall, mapvar, vmin=v_min, vmax=v_max, cmap=color_code)  #plt.cm.hot_r)
  #themapc = plt.contour(lonall, latall, mapvar, 10, color=(0,0,0))
  #plt.clabel(themapc, fontsize=12)
  #themap.cmap.set_over((0.0,0.0,0.0))
  #themap.cmap.set_under((0.5,0.5,0.5))
  plt.title(title, fontsize=12)
  plt.axis([0, 356.5, -90, 90])
  #plt.xlabel('Longitude [' + lon_units + ']')
  #plt.ylabel('Latitude [' + lat_units + ']')
  plt.colorbar(themap, orientation='horizontal', extend='both', ticks=[-200,-100,-50,-20,-10,-5,-2,-1,-0.5,0,0.5,1,2,5,10,20,50,100,200])
  plt.savefig(outfile)
  plt.draw()
  if Show:
    plt.show()
  return

