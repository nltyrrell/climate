import numpy as np
import Scientific.IO.NetCDF as S

def algorhythm_add(nc_anom,nc_clim='ncfiles/check.cxf.sst.nc',nc_out='testy.nc',scale=1):
#  txt_file='int_corr.txt'
#  anomaly = readtext(txt_file)
#  print 'Text file read in'
  
#  nc_file='ncfiles/reg0lag.nc'
  variable='temp'
  anomaly = readnetcdf(nc_anom,variable)
  print 'Netcdf file read in'

#  nc_file='ncfiles/check.cxf.sst.nc'
  variable='temp'
  clim = readnetcdf(nc_clim,variable)
  print 'Netcdf file read in'
  
  ancil = add_anom(clim,anomaly,scale)
  print 'Ancil array created'
  print 'Scale =',scale
  
  writenc(nc_out,ancil)
  print 'Netcdf written'
  return anomaly, ancil

def readtext(txt_file):
  # Function reads in correctly formatted text file of SST anomalies, outputs numpy array
  fileobj = open(txt_file, 'r')
  outputstr = fileobj.readlines()
  fileobj.close()
  textarray = np.array([[float(n) for n in line.split()] for line in outputstr])
  return textarray

def readnetcdf(nc_file,variable):
  """ Function the reads in an nc file and outputs 
  a numpy array (I think)
  """
  nc_fileobj = S.NetCDFFile(nc_file,mode='r')
  nc_data = nc_fileobj.variables[variable].getValue()[:,0,:,:]
  return nc_data

def add_anom(climatology, anomaly, scale):
  # Adds an anomaly to a climatology
  #repdim_anom = np.repeat(np.expand_dims(anomaly,axis=0),12,0)
  ancil=np.zeros([12,73,96])
  ancil=climatology+anomaly*scale
  return ancil

def writenc(nc_fileout,ancil):
  nc_out = S.NetCDFFile(nc_fileout,mode='a')
  var=nc_out.variables['temp']
  data_out=var.getValue()
  ancil2=np.expand_dims(ancil,axis=1)
  data_out[:,:,:,:]=ancil2[:,:,:,:]
  var[:]=data_out
  nc_out.close()
  return ancil2


#  print(np.shape(climatology))
#  print(np.shape(anomaly))
#  for i in xrange(12):
#    ancil=np.zeros([12,73,96])
#    print(i)
#    ancil[i,:,:] = climatology[i,:,:] + anomaly[:,:]  #*scale
