import numpy as np
import inout as io

def troposave(data):
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


	meanT = (data[:,3:,::]+data[:,2:-1,::])*0.5

	meanRho = (density[3:]+density[2:-1])*0.5

	diffH = np.diff(height[2:])
	massw = meanRho * diffH
	w_ave = np.average(meanT,axis=1,weights=massw)

	return w_ave


