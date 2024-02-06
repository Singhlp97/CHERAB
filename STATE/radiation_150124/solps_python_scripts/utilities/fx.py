import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

from solps_python_scripts.read_b2fgmtry import read_b2fgmtry

def compute_bp(ix = None, iy = None, bb = None):

	return bb[ix,iy:,0] # np.sqrt(bb[ix,iy:,0]**2 + bb[ix,iy:,1]**2)

def compute_r(ix = None, iy = None, crx = None):

	rx = []
	for y in range(iy, crx.shape[1]): rx.append(crx[ix,y,:].mean())

	return np.array(rx)

def compute_fx(where = ".", save = True):

	b2fgmtry = read_b2fgmtry(where = where)

	nx = b2fgmtry['nx']
	ny = b2fgmtry['ny']

	bb = b2fgmtry['bb']
	crx = b2fgmtry['crx']

	# CAUTION.
	#
	# - first  column = inner target values
	# - second column = outer target values
	# - double null   = two targets only 
	# - mid-planes mid-way through left-/right-cuts
	#
	# TODO:
	#
	# - jxi and jxa from b2mn.dat instead of left-/right-cuts

	# - only SOL region need
	# - assuming topcut[0] == topcut[1]
	iy = int(b2fgmtry['topcut'][0])

	fx    = np.zeros((int(iy), 2))
	bpmp  = np.zeros((int(iy), 2))
	bpdiv = np.zeros((int(iy), 2))
	rmp   = np.zeros((int(iy), 2))
	rdiv  = np.zeros((int(iy), 2))

	########

	# inner bottom target: ix = 0
	ix = 0
	bpdiv[:,0] = compute_bp(ix, iy, bb)
	rdiv[:,0]  = compute_r(ix, iy, crx)

	# inner midplane: ix = mid-way through left cuts
	ix = int(5E-01 * b2fgmtry['leftcut'].sum())
	bpmp[:,0] = compute_bp(ix, iy, bb)
	rmp[:,0]  = compute_r(ix, iy, crx)

	# outer midplane: ix = mid-way through right cuts
	ix = int(5E-01 * b2fgmtry['rightcut'].sum())
	bpmp[:,1] = compute_bp(ix, iy, bb)
	rmp[:,1]  = compute_r(ix, iy, crx)

	# outer top target: ix = -1
	ix = -1
	bpdiv[:,1] = compute_bp(ix, iy, bb)
	rdiv[:,1]  = compute_r(ix, iy, crx)

	########

	# fx computation

	fx[:,0] = rmp[:,0] * bpmp[:,0] / rdiv[:,0] / bpdiv[:,0]
	fx[:,1] = rmp[:,1] * bpmp[:,1] / rdiv[:,1] / bpdiv[:,1]

	########

	# save data

	if save is True:

		data = {}

		data['rmp']   = rmp
		data['rdiv']  = rdiv
		data['bpmp']  = bpmp
		data['bpdiv'] = bpdiv
		data['fx']    = fx

		pickle.dump(data, open(os.path.join(where, 'fx_data.pkl'), 'wb'))

	return data