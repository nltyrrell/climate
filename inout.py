import numpy as np
import scipy as sp
import Scientific.IO.NetCDF as S

# Import netcdf file as array
def readnc(nc_file,variable):
  """ Function the reads in an nc file and outputs a numpy array """
  nc_fileobj = S.NetCDFFile(nc_file,mode='r')
  nc_data = nc_fileobj.variables[variable].getValue()
  return nc_data, nc_fileobj

# Export result map
def writenc4d(nc_fileout,ancil):
  """ Use regmap, or similar, and add anomaly to a 12 month climatology """
  nc_out = S.NetCDFFile(nc_fileout,mode='a')
  var=nc_out.variables['temp']
  data_out=var.getValue()
  #anc_exp = np.repeat(np.expand_dims(ancil,axis=0),12,0)
  #ancil2=np.expand_dims(anc_exp,axis=1)
  data_out[:,:,:,:]=ancil[:,:,:,:]
  var[:]=data_out
  nc_out.close()
  return ancil

# Import timeseries
def readtext(txt_file,col=0):
	# Function reads in correctly formatted text file of SST anomalies, outputs numpy array
	fileobj = open(txt_file, 'r')
	outputstr = fileobj.readlines()
	fileobj.close()
	textarray = np.array([[float(n) for n in line.split()] for line in outputstr])[:,col]
	return textarray

def writetext(data,text_file):
	#could just use np.savetext('textfile.txt',data,fmt='%10.5f')
	fileout = open(text_file, 'w') 
	fileout.writelines(data) 
	fileout.close() 



# Export result map
def writenc(nc_fileout,outdata,var,dims=4,ext_data=None):
	""" input: nc_fileout - nc file that will be overwritten
	outdata - a (t,z,lat,lon) array to write to ncfile
	var - e.g. temp, sst of ncfile
	dims - dimensions of outdata
	"""
	nc_out = S.NetCDFFile(nc_fileout,mode='a')
	values=nc_out.variables[var]
	data_out=values.getValue()
	#if dims==4:
	#  ancil_out=np.expand_dims(ancil,axis=1)
	#elif dims==3:
	if ext_data!=None:
		extd = np.zeros(data_out.shape)
		outdata = extd + outdata
	data_out[:,:,:,:]=outdata[:,:,:,:]
	values[:]=data_out
	nc_out.close()
	return 





