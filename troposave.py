import numpy as np
import iris
 
def troposave(data, startlev=0, endlev=11,weights=False):
    # pressure in Pa
    plev = np.array([100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000,\
                    20000, 15000, 10000, 7000, 5000, 3000, 2000])
    # geometr. height in m
    height = np.array([110.9, 762.1, 1457.6, 3013.6, 4209.2, 5579.3, 7193.6, 9177.2, 10379.9,\
            11805.9, 13637.6, 16221, 18495.3, 20643, 23938.4, 26592.3])
    # temperature in C
    temperature = np.array([14.28, 10.05, 5.53, -4.58, -12.34, -21.23, -31.71, -44.57, -52.36,\
            -56.5, -56.5, -56.5, -56.5, -55.92, -52.65, -50.02])
    # density in kg/cm
    density = np.array([1.212, 1.1378, 1.0625, 0.90796, 0.80142, 0.69142, 0.57713, 0.4572, 0.39445,\
            0.32159, 0.24119, 0.16079, 0.11256, 8.0184e-2, 4.7397e-2, 3.1223e-2])
    # Calulate a mass weighted tropospheric mean


    meanT = (data[:,startlev+1:endlev,::]+data[:,startlev:endlev-1,::])*0.5

    meanRho = (density[startlev+1:endlev]+density[startlev:endlev-1])*0.5

    diffH = np.diff(height[startlev:endlev])
    massw = meanRho * diffH
    w_ave = np.average(meanT,axis=1,weights=massw)
    if weights == True:
        outdata = [w_ave, massw]
    else:
        outdata = w_ave
    return outdata

def troposweights(startlev=0, endlev=11):
    # pressure in Pa
    plev = np.array([100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000,\
                    20000, 15000, 10000, 7000, 5000, 3000, 2000])
    # geometr. height in m
    height = np.array([110.9, 762.1, 1457.6, 3013.6, 4209.2, 5579.3, 7193.6, 9177.2, 10379.9,\
            11805.9, 13637.6, 16221, 18495.3, 20643, 23938.4, 26592.3])
    # temperature in C
    temperature = np.array([14.28, 10.05, 5.53, -4.58, -12.34, -21.23, -31.71, -44.57, -52.36,\
            -56.5, -56.5, -56.5, -56.5, -55.92, -52.65, -50.02])
    # density in kg/cm
    density = np.array([1.212, 1.1378, 1.0625, 0.90796, 0.80142, 0.69142, 0.57713, 0.4572, 0.39445,\
            0.32159, 0.24119, 0.16079, 0.11256, 8.0184e-2, 4.7397e-2, 3.1223e-2])
    # Calulate a mass weighted tropospheric mean


    #meanT = (data[:,startlev+1:endlev,::]+data[:,startlev:endlev-1,::])*0.5

    meanRho = (density[startlev+1:endlev]+density[startlev:endlev-1])*0.5

    diffH = np.diff(height[startlev:endlev])
    massw = meanRho * diffH
    #w_ave = np.average(meanT,axis=1,weights=massw)

    return massw

def tropwtsiris(cube):
    """ 
    Calculate the mass weights for a troposhperic average
    Input:  Iris cube, constrained in pressure levels must use 'air_pressure'
    Output: broadcast weights, mass weights
    """

    # pressure in Pa
    pressure = np.array([ 1000.,   925.,   850.,   700.,   600.,   500.,   400.,   300.,
                        250.,   200.,   150.,   100.,    70.,    50.,    30.,    20.]) 
    plev = cube.coord('air_pressure').points

    print plev[0]
    print plev[-1]
    startlev = np.where(pressure == np.round(plev[0]))[0][0]
    endlev = np.where(pressure == np.round(plev[-1]))[0][0] + 2
    print startlev
    print endlev

    # geometr. height in m
    height = np.array([110.9, 762.1, 1457.6, 3013.6, 4209.2, 5579.3, 7193.6, 9177.2, 10379.9,\
            11805.9, 13637.6, 16221, 18495.3, 20643, 23938.4, 26592.3])
    # temperature in C
    temperature = np.array([14.28, 10.05, 5.53, -4.58, -12.34, -21.23, -31.71, -44.57, -52.36,\
            -56.5, -56.5, -56.5, -56.5, -55.92, -52.65, -50.02])
    # density in kg/cm
    density = np.array([1.212, 1.1378, 1.0625, 0.90796, 0.80142, 0.69142, 0.57713, 0.4572, 0.39445,\
            0.32159, 0.24119, 0.16079, 0.11256, 8.0184e-2, 4.7397e-2, 3.1223e-2])
    # Calulate a mass weighted tropospheric mean


    #meanT = (data[:,startlev+1:endlev,::]+data[:,startlev:endlev-1,::])*0.5

    meanRho = (density[startlev+1:endlev]+density[startlev:endlev-1])*0.5

    diffH = np.diff(height[startlev:endlev])
    massw = meanRho * diffH
    print cube.shape
    print massw.shape
    #w_ave = np.average(meanT,axis=1,weights=massw)
    bcweight = iris.util.broadcast_to_shape(massw,cube.shape,(1,))
    return bcweight, massw

def iristropave(cube, plev_bottom=1000, plev_top=200):
    """
    Calculate mass weighted tropospheric average
    Input:  Iris cube, bottom pressure level, top pressure level
    Output: Iris cube, one vertical level
    Units: hPa, 'air_pressure'
    """
    
    trop = iris.Constraint(air_pressure = lambda p: plev_top <= p <= plev_bottom)
    cube_trop = cube.extract(trop)

    bcweight = tropwtsiris(cube_trop)[0]

    trop_weighted = cube_trop.collapsed('air_pressure',
                               iris.analysis.MEAN,
                               weights=bcweight)

    return trop_weighted






