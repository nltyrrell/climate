import Scientific.IO.NetCDF as S
import numpy as np
import matplotlib.pyplot as plt
import sys
import mpl_toolkits.basemap as bm
import scipy as sc
from scipy import stats

# Read in the netcdf files

def readnetcdf(nc_file,variable):
  """ Function the reads in an nc file and outputs 
  a numpy array (I think, fairly sure)
  Input: nc_data = readnetcdf('ncfiles/file.nc','temp')
  """
  nc_fileobj = S.NetCDFFile(nc_file,mode='r')
  nc_data = nc_fileobj.variables[variable].getValue()[:,0,:,:]
  return nc_data, nc_fileobj

#clim = S.NetCDFFile('~/project/access/clim/ncfiles/temp.sfc.clim.nc',mode='r')

def mean_ttest(nc_data,nc_clim,sig=0.05):
  """ Calculate the mean, minus climatology (i.e. get anomalies), 
  use ttest and mask insignificant areas 
  """
  var_mean  = nc_data.mean(axis=0)
  clim_mean  = nc_clim.mean(axis=0)
  var_anom = var_mean[:,:] - clim_mean[:,:]
  # Calculate the t-test, get rid of useless height axis
  tt_var = sc.stats.ttest_rel(nc_clim[0:150,:,:],nc_data[0:150,:,:],axis=0)
  # Mask out insignificant areas
  tt_masked = np.ma.masked_where(tt_var[1]>=0.05,var_anom)
  return var_anom, tt_var, tt_masked

def makemap(nc_obj,mapvar,color_code,v_min=None,v_max=None):
  """
  Map a climate variable using basemap and matplotlib
  for color_code, need to do plt.cm.hot_r etc
  """
  
  if v_min == None:
    v_min=mapvar.mean()-2*mapvar.std()
    print v_min
  elif v_max == None:
    v_max=mapvar.mean()+2*mapvar.std()
    print v_max

  lon=nc_obj.variables['longitude'].getValue()
  lon_units = nc_obj.variables['longitude'].units
  lat = nc_obj.variables['latitude'].getValue()
  lat_units = nc_obj.variables['latitude'].units
  #- Make 2-D longitude and latitude arrays:
  [lon2d, lat2d] = np.meshgrid(lon, lat)

  #- Set up map:
  plt.clf()
  mapproj = bm.Basemap(projection='cyl', llcrnrlat=-90.0, llcrnrlon=0.0, urcrnrlat=90.0, urcrnrlon=360.0)
  mapproj.drawcoastlines()
  mapproj.drawparallels(np.array([-90, -45, 0, 45, 90]), labels=[1,0,0,0])
  mapproj.drawmeridians(np.array([0, 90, 180, 270, 360]), labels=[0,0,0,1])
  lonall, latall = mapproj(lon2d, lat2d)

  themap = plt.pcolormesh(lonall, latall, mapvar, vmin=v_min, vmax=v_max, cmap=color_code)  #plt.cm.hot_r)
  themapc = plt.contour(lonall, latall, mapvar, 10, color=(0,0,0))
  plt.clabel(themapc, fontsize=12)
  themap.cmap.set_over((0.0,0.0,0.0))
  themap.cmap.set_under((1.0,1.0,1.0))
  plt.title('DJF Mean Temperature response, [K]', fontsize=12)
  plt.axis([0, 360, -90, 90])
  #plt.xlabel('Longitude [' + lon_units + ']')
  #plt.ylabel('Latitude [' + lat_units + ']')
  plt.colorbar(themap, orientation='horizontal', extend='both', ticks=[-200,-100,-50,-20,-10,-5,-2,-1,0,1,2,5,10,20,50,100,200])
  plt.savefig('temp.test.eps')
  plt.show()
  return
  #raise SystemExit


def algorhythm(ncfile,ncclim,variable, vmin, vmax):
  
  clim_data,clim_obj = readnetcdf(ncclim,variable)
  nc_data,nc_obj = readnetcdf(ncfile,variable)
  var_anom,tt_var,tt_masked  = mean_ttest(nc_data,clim_data)

  makemap(nc_obj, tt_masked, plt.cm.hot_r, vmin, vmax)

  return





