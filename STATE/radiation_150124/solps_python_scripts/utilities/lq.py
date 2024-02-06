import warnings

import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

from solps_python_scripts.read_b2mn import read_b2mn
from solps_python_scripts.utilities.last10s import read_last10s

############################################################################################

def exp(r, q0, lq):
	return q0 * np.exp(- r / lq)

############################################################################################

def compute_lq(where = ".", num_points = None, legend = False):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	if not isinstance(where, list): where = [where]

	for here in where:

		last10s = read_last10s(where = here, save = True)

		# lq as computed by remapping ft3dr at on dsa:
		# - ft3dr (and dsa) from peak onwards

		b2mn = read_b2mn(where = here)

		ft = last10s['ft3drP'][:,1]
		dsa = last10s['dsa']
		start = np.where(ft == ft.max())[0][0]

		ft = ft[start:]
		if int(b2mn['b2mwti_jxi']) == 0: dsa = dsa[dsa > 0]
		dsa = dsa[start:]

		if num_points is not None:
			dsa = dsa[:num_points]
			ft  = ft[:num_points]

		param0 = (ft.max(), 5E-03)
		param, cv = scipy.optimize.curve_fit(exp, dsa, ft, param0)
		q0, lq = param

		squaredDiffs = np.square(ft - exp(dsa, q0, lq))
		squaredDiffsFromMean = np.square(ft - np.mean(ft))
		rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)

		plt.figure()
		plt.plot(dsa, ft, 'ko-', label = here)
		plt.plot(dsa, exp(dsa, q0, lq), 'm--')
		plt.xlabel('$s\;[m]$')
		plt.ylabel('$[MW \cdot m^{-2}]$')
		plt.title('$\lambda_q =$ {:.3G} mm - $R^2 =$ {:.3G}'.format(lq * 1E+03, rSquared))
		plt.grid(True)
		if legend is True: plt.legend(loc = 'best')

	plt.show()

	if rSquared < 9.5E-01:
		warnings.warn('R^2 of exponential fit is of poor quality! Check lq.py!')

	return 

############################################################################################

def compute_ln(where = ".", num_points = None, legend = False):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	if not isinstance(where, list): where = [where]

	for here in where:

		last10s = read_last10s(where = here, save = True)

		# lq as computed by remapping ft3dr at on dsa:
		# - ft3dr (and dsa) from peak onwards

		b2mn = read_b2mn(where = here)

		keys = ['fl3drP', 'fo3drP']
		(dsa, fn, q0, ln, rSquared)	= compute_fit(keys = keys,
												  last10s = last10s,
												  b2mn = b2mn,
												  num_points = num_points)

		for i in range(dsa.shape[1]):

			plt.figure()
			plt.plot(dsa[:,i], fn[:,i], 'ko-', label = here, linewidth = 3)
			plt.plot(dsa[:,i], exp(dsa[:,i], q0[i], ln[i]), 'm-', linewidth = 1)
			plt.xlabel('$s\;[m]$')
			plt.ylabel('$[m^{-2}]$')
			if keys[i] == 'fl3drP': species = 'electrons'
			if keys[i] == 'fo3drP': species = 'ions'
			plt.title('{}: $\lambda_n =$ {:.3G} mm - $R^2 =$ {:.3G}'.format(species, ln[i] * 1E+03, rSquared[i]))
			plt.grid(True)
			if legend is True: plt.legend(loc = 'best')

		plt.show()

############################################################################################

def compute_fit(keys = None, last10s = None, b2mn = None, num_points = None):

	ff = []
	ds = []
	ll = []
	qq0 = []
	rSquared = []

	for (i,key) in enumerate(keys):

		f = last10s[key][:,1]
		dsa = last10s['dsa']
		start = np.where(f == f.max())[0][0]

		f = f[start:]
		if int(b2mn['b2mwti_jxi']) == 0: dsa = dsa[dsa > 0]
		dsa = dsa[start:]

		if num_points is not None:
			dsa = dsa[:num_points]
			f   = f[:num_points]

		param0 = (f.max(), 5E-03)
		param, cv = scipy.optimize.curve_fit(exp, dsa, f, param0)
		q0, l = param

		squaredDiffs = np.square(f - exp(dsa, q0, l))
		squaredDiffsFromMean = np.square(f - np.mean(f))

		ff += [f]
		ds += [dsa]
		ll += [l]
		qq0 += [q0]
		rSquared += [1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)]

	if (np.array(rSquared) < 9.5E-01).sum() > 0:
		warnings.warn('R^2 of exponential fit is of poor quality! Check lq.py!')

	return (np.array(ds).T, np.array(ff).T, np.array(qq0), np.array(ll), np.array(rSquared))