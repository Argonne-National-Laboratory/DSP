from scipy.interpolate import interp1d
import numpy as np

for n in [1, 6, 12, 48, 96]:
	for k in [1,2,3]:

		f = open('wp64-' + str(k) + '.txt')

		m = np.matrix([[float(i) for i in line.split()] for line in f.readlines()])

		cols = []

		for j in range(m.shape[1]):
			y = np.transpose(m[:,j])
			intf = interp1d(np.arange(y.shape[1]), y, kind='cubic')
			new_col = np.array([intf(i)[0] for i in np.linspace(0, y.shape[1]-1, n)])
			cols.append(new_col)

		f.close()

		fo = open('wp64-' + str(n) + '-' + str(k) + '.txt', 'w')
		for j in range(n):
			for i in range(len(cols)):
				fo.write(str(cols[i][j]) + ' ')
			fo.write('\n')
		fo.close()