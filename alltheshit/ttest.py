import Scientific.IO.NetCDF as S
import numpy as np
import matplotlib.pyplot as plt
import sys
import mpl_toolkits.basemap as bm
import scipy as sc
from scipy import stats
from mycmaps import jetwhite

# Read in the netcdf files

def readnetcdf(nc_file,variable):
  """ Function the reads in an nc file and outputs 
  a numpy array
  Input: nc_data = readnetcdf('ncfiles/file.nc','temp')
  """
  nc_fileobj = S.NetCDFFile(nc_file,mode='r')
  nc_data = nc_fileobj.variables[variable].getValue()[:,0,:,:]
  return nc_data, nc_fileobj

def mean_ttest(nc_data,nc_clim,sig=0.05):
  """ Calculate the mean, minus climatology (i.e. get anomalies)
  use ttest and mask insignificant areas 
  """
  var_mean  = nc_data.mean(axis=0)
  clim_mean  = nc_clim.mean(axis=0)
  var_anom = var_mean[:,:] - clim_mean[:,:]
  # Calculate the t-test, get rid of useless height axis
  tt_var = sc.stats.ttest_rel(nc_clim[0:360,:,:],nc_data[0:360,:,:],axis=0)
  # Mask out insignificant areas
  tt_masked = np.ma.masked_where(tt_var[1]>=0.05,var_anom)
  return var_anom, tt_var, tt_masked

def makemap(nc_obj,mapvar,color_code=jetwhite(n=256),vmin=None,vmax=None, name='test', title='title'):
  """
  Map a climate variable using basemap and matplotlib
  for color_code, need to do plt.cm.hot_r etc
  """
  
  if vmin == None:
    vmin=np.floor(mapvar.mean()-2*mapvar.std())
    print vmin
  elif vmax == None:
    vmax=np.floor(mapvar.mean()+2*mapvar.std())
    print vmax

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

  plt.ion()
  themap = plt.pcolormesh(lonall, latall, mapvar, vmin=vmin, vmax=vmax, cmap=color_code)  #plt.cm.hot_r)
  #themap = plt.contourf(lonall, latall, mapvar, np.linspace(vmin,vmax,21), cmap=color_code, extend="both")
  #plt.clabel(themapc, fontsize=12)
  themap.cmap.set_over((0.0,0.0,0.0))
  themap.cmap.set_under((1.0,1.0,1.0))
  plt.title(title, fontsize=14)
  plt.axis([0, 357, -90, 90])
  #plt.xlabel('Longitude [' + lon_units + ']')
  #plt.ylabel('Latitude [' + lat_units + ']')
  plt.colorbar(themap, orientation='horizontal', extend='both', ticks=np.linspace(vmin,vmax,11))
  plt.savefig(name + '.eps')
  plt.show()
  return
  #raise SystemExit


def algorhythm(ncfile,ncclim,variable, vmin=-1, vmax=1, title='test',name='test'):
  
  clim_data,clim_obj = readnetcdf(ncclim,variable)
  nc_data,nc_obj = readnetcdf(ncfile,variable)
  var_anom,tt_var,tt_masked  = mean_ttest(nc_data,clim_data)

  makemap(nc_obj, var_anom, vmin=vmin, vmax=vmax,title=title,name=name)

  return var_anom, nc_obj, 

#if __name__ = '__main__':
exper=['r0p1','r6g2','r6p2','r6r2']
#title_list=['Pressure anomalies, Pacific 0 lag pattern, [Pa]','Pressure anomalies, Tropical 6m lag pattern, [Pa]','Pressure anomalies, Pacific 6m lag pattern, [Pa]','Pressure anomalies, Atl/Ind 6m lag pattern, [Pa]']
#title_list=['Longwave down anomalies, Pacific 0 lag pattern, [Wm-2]','Longwave down anomalies, Tropical 6m lag pattern, [Wm-2]','Longwave down anomalies, Pacific 6m lag pattern, [Wm-2]','Longwave down anomalies, Atl/Ind 6m lag pattern, [Wm-2]']
n=0
outvar=[1,2,3,4]
var=['temp','pres','dlwf','dswf']
varname=['temp','p','ilr','field208']
#Choose variable
v=0
for i in exper:
  ncfile='ncfiles/'+ var[v] +'.sfc.'+ i +'.nc'
  ncclim='../clim/ncfiles/'+ var[v] +'.sfc.clim.nc'
  variable=varname[v]
  title=var[v] +' '+ exper[n] #title_list[n]
  name=var[v]+'.'+exper[n]
  outvar[n] = algorhythm(ncfile,ncclim,variable,vmin=-1,vmax=1,title=title,name=name)

  n=n+1




