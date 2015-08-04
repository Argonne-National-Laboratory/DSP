from scipy.interpolate import interp1d
import numpy as np

y = [2389.12, 2239.80, 2165.14, 2090.48, 2090.48, 2165.14,
	 2389.12, 2837.08, 3247.71, 3546.35, 3695.67, 3733.00,
	 3695.67, 3733.00, 3733.00, 3621.01, 3583.68, 3583.68,
	 3471.69, 3434.36, 3434.36, 3471.69, 3247.71, 2687.76]
f = interp1d(np.arange(len(y)), y, kind='cubic')

n = 96

for i in np.linspace(0, len(y)-1, n):
	print str(f(i)) + '\t0.03\t0.03'
