from cdo import *
import numpy as np
import scipy as sp
import shutil as sh

cdo=Cdo()

def repeat(infile,repeats,nc_length=1)
  """ enter the file to be repeated, eg clim, how many repeats (one less year than total required length), optional: length of the file, normally 1 year, unit of time to repeat over, normally days.
  example: repeat(infile='clim.nc',repeats='25',nc_length='4') """

  sh.copyfile(infile,'infile.nc')
  
  for i in np.arange(24):
    print(i+1)
    timeshift=nc_length*(i+1)
    cdo.shifttime(str(timeshift)+'years',output='shifted.nc',input=infile,options='-O')
    outfile=cdo.mergetime(output='outfile.nc',input='infile.nc shifted.nc',options='-O')
    sh.copyfile('outfile.nc','infile.nc')
    print(i+1)

  return




