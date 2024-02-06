import os
import time
import pickle
import numpy as np

from solps_python_scripts.read_b2mn import read_b2mn
from solps_python_scripts.read_b2fgmtry import read_b2fgmtry
from solps_python_scripts.read_b2fstate import read_b2fstate
from solps_python_scripts.read_b2fplasmf import read_b2fplasmf
from solps_python_scripts.read_ft44 import read_ft44

from solps_python_scripts.reactions.compute_rates import compute_rates

from solps_python_scripts.utilities.accessories import load_pickle

#########################################################################################################

def read_last10s(where = ".", extended = False, save = True):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	try:    time_pkl = os.path.getmtime(os.path.join(where, 'last10s.pkl'))
	except: time_pkl = 0

	try:    time_runlog = os.path.getmtime(os.path.join(where, 'run.log'))
	except: time_runlog = 1

	if time_runlog < time_pkl:

		print('last10s loaded from .pkl')
		return pickle.load(open(os.path.join(where, 'last10s.pkl'), 'rb'))

	else:

		print('read_last10s')

		original = os.getcwd()
		os.chdir(where)

		if extended is True: os.system("2d_profiles_extended &")
		else:                os.system("2d_profiles &")
		time.sleep(3.0)

		b2mn = read_b2mn(where = '.')

		last10s = {}

		for last10 in os.listdir():
			ext = last10.split(".")[-1]
			if ext == "last10" and last10.split(".")[0] != "run":
				if os.stat(last10).st_size > 0:
					# try:
					if int(b2mn['b2mwti_jxi']) == 0 and last10.split(".")[0][-1] == 'i':
						# limiter => different convention
						data = np.loadtxt(open(last10, "r"))
						last10s[last10.split(".")[0]] = data[np.where(data[:,0] < 0)]
					else:
						# usual convention
						last10s[last10.split(".")[0]] = np.loadtxt(open(last10, "r"))
					# except: pass
		
		os.chdir(original)

		# last10s['fe3dl'][:,1] *= 1E-06
		# last10s['fe3dr'][:,1] *= 1E-06
		# last10s['fi3dl'][:,1] *= 1E-06
		# last10s['fi3dr'][:,1] *= 1E-06
		# last10s['ft3dl'][:,1] *= 1E-06
		# last10s['ft3dr'][:,1] *= 1E-06
		# last10s['fe3dlP'][:,1] *= 1E-06
		# last10s['fe3drP'][:,1] *= 1E-06
		# last10s['fi3dlP'][:,1] *= 1E-06
		# last10s['fi3drP'][:,1] *= 1E-06
		# last10s['ft3dlP'][:,1] *= 1E-06
		# last10s['ft3drP'][:,1] *= 1E-06
		# try:
		# 	last10s['fe3dtl'][:,1] *= 1E-06
		# 	last10s['fe3dtr'][:,1] *= 1E-06
		# 	last10s['fi3dtl'][:,1] *= 1E-06
		# 	last10s['fi3dtr'][:,1] *= 1E-06
		# 	last10s['ft3dtl'][:,1] *= 1E-06
		# 	last10s['ft3dtr'][:,1] *= 1E-06
		# except:
		# 	pass

		# add further data

		last10s = add_ds(last10s = last10s, where = where)
		last10s = add_vol(last10s = last10s, where = where)
		try: last10s = add_sheath(last10s = last10s)
		except: pass
		# last10s = add_sheath(last10s = last10s)
		if b2mn['b2mndr_eirene'] == 1: last10s = add_ft44_b2(last10s = last10s, where = where)
		last10s = add_na(last10s = last10s, where = where)
		last10s = add_Zeff(last10s = last10s, where = where)
		last10s = add_transport(last10s = last10s, where = where)
		if b2mn['b2mndr_eirene'] == 1: last10s = add_rates(last10s = last10s, where = where)
		# last10s = add_grazing_angle(last10s = last10s, where = where)
		# last10s = add_fht(last10s = last10s, where = where)

		if save is True:
			pickle.dump(last10s, open(os.path.join(where, 'last10s.pkl'), 'wb'))

		return last10s

#########################################################################################################

def add_ds(where = ".", last10s = None):

	for location in ["a", "i", "l", "r"]:
		try: last10s["ds" + location] = np.loadtxt(open(os.path.join(where, "ds" + location), "r"))
		except: pass

	return last10s

#########################################################################################################

def add_vol(where = ".", last10s = None):

	# Author: Matteo Moscheni
	# E-mail: matteo.moscheni@tokamakenergy.co.uk
	# February 2022

	b2mn = read_b2mn(where = where)
	b2fgmtry = read_b2fgmtry(where = where, verbose = False)

	locations = ['a', 'i', 'l', 'r']

	for location in locations:

		data = b2fgmtry['vol']

		label = "".join(['vol', '3d', location])

		last10 = create_profile(b2mn = b2mn, b2fgmtry = b2fgmtry,
								last10s = last10s, data = data, location = location)

		last10s[label] = last10

	return last10s

#########################################################################################################

def create_profile(b2mn = None, b2fgmtry = None, last10s = None, data = None, location = None):

	if location == 'a':
		# ix = int(b2mn['b2mwti_jxa'])
		if b2fgmtry['nncut'][0] == 0:
			ix = int(b2fgmtry['nx'] / 2)
		elif b2fgmtry['nncut'][0] == 1:
			ix = int(b2fgmtry['rightcut'][0] - (b2fgmtry['rightcut'][0] - b2fgmtry['leftcut'][0]) / 4)
		elif b2fgmtry['nncut'][0] == 2:
			ix = int((b2fgmtry['rightcut'][0] + b2fgmtry['rightcut'][1]) / 2)
		# first ghost cell in ix=0 has been deleted
		# and start from 0 not from -1
		ix += 2
	elif location == 'i':
		# ix = int(b2mn['b2mwti_jxi'])
		if b2fgmtry['nncut'][0] == 0:
			ix = 0
		elif b2fgmtry['nncut'][0] == 1:
			ix = int(b2fgmtry['leftcut'][0] + (b2fgmtry['rightcut'][0] - b2fgmtry['leftcut'][0]) / 4)
		elif b2fgmtry['nncut'][0] == 2:
			ix = int((b2fgmtry['leftcut'][0] + b2fgmtry['leftcut'][1]) / 2)
		# first ghost cell in ix=0 has been deleted
		# and start from 0 not from -1
		ix += 2
	elif location == 'l':
		ix = 0
	elif location == 'r':
		ix = -1

	try:
		n = data.shape[2]
		datum = data[ix,:,:]
	except:
		n = 1
		datum = data[ix,:]

	try:    abscissa = last10s["".join(['ds', location])][1:-1]
	except: abscissa = last10s["".join(['ne3d', location])][1:-1]

	if int(b2mn['b2mwti_jxi']) == 0 and location != 'a':
		# limiter case => different convention
		if location == 'i':
			if n == 1:
				last10 = np.column_stack((abscissa, datum[:len(abscissa)]))
			else:
				last10 = np.column_stack((abscissa, datum[:len(abscissa),:]))
		elif location == 'l' or location == 'r':
			if n == 1:
				last10 = np.column_stack((abscissa, datum[-len(abscissa):]))
			else:
				last10 = np.column_stack((abscissa, datum[-len(abscissa):,:]))
	else:
		# single null and double null
		last10 = np.column_stack((abscissa, datum))

	return last10


#########################################################################################################

def add_sheath(last10s = None):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	ti3dl = last10s['ti3dl'][:,1]
	ti3dr = last10s['ti3dr'][:,1]
	te3dl = last10s['te3dl'][:,1]
	te3dr = last10s['te3dr'][:,1]

	po3dl = last10s['po3dl'][:,1]
	po3dr = last10s['po3dr'][:,1]

	l = last10s['ti3dl'][:,0]
	r = last10s['ti3dr'][:,0]

	# default values

	gi0 = 1.5 * np.ones(last10s['ti3dl'].shape[0])
	ge0 = 1.0 * np.ones(last10s['ti3dl'].shape[0])

	# internal energy

	last10s['sei3dl'] = np.array([l, ge0 + po3dl / te3dl]).T
	last10s['sei3dr'] = np.array([r, ge0 + po3dr / te3dr]).T
	last10s['sii3dl'] = np.array([l, gi0]).T
	last10s['sii3dr'] = np.array([l, gi0]).T

	# tei = Ti/Te

	last10s['tie3da'] = np.array([last10s['ti3da'][:,0], last10s['ti3da'][:,1]/last10s['te3da'][:,1]]).T
	last10s['tie3di'] = np.array([last10s['ti3di'][:,0], last10s['ti3di'][:,1]/last10s['te3di'][:,1]]).T
	last10s['tie3dl'] = np.array([last10s['ti3dl'][:,0], last10s['ti3dl'][:,1]/last10s['te3dl'][:,1]]).T
	last10s['tie3dr'] = np.array([last10s['ti3dr'][:,0], last10s['ti3dr'][:,1]/last10s['te3dr'][:,1]]).T

	# total energy
	# - inertialess electrons
	# - default cs expression

	last10s['set3dl'] = last10s['sei3dl']
	last10s['set3dr'] = last10s['sei3dl']
	last10s['sit3dl'] = np.array([l, gi0 + 5E-01 * (1.0 + te3dl / ti3dl)]).T
	last10s['sit3dr'] = np.array([l, gi0 + 5E-01 * (1.0 + te3dr / ti3dr)]).T

	last10s['stt3dl'] = np.array([l, last10s['set3dl'][:,1] + last10s['sit3dl'][:,1]]).T
	last10s['stt3dr'] = np.array([l, last10s['set3dr'][:,1] + last10s['sit3dr'][:,1]]).T

	# simil-kinetic

	# ion internatml energy transmission

	last10s['siik3dl'] = np.array([l, 5E-01 * (5 * ti3dl / te3dl - te3dl / ti3dl - 1.0)]).T
	last10s['siik3dr'] = np.array([l, 5E-01 * (5 * ti3dr / te3dr - te3dr / ti3dr - 1.0)]).T

	# ion total energy transmission

	last10s['sitk3dl'] = np.array([l, 2.5E-00 * ti3dl / te3dl]).T
	last10s['sitk3dr'] = np.array([l, 2.5E-00 * ti3dr / te3dr]).T

	# total energy transmission

	last10s['sttk3dl'] = np.array([l, last10s['set3dl'][:,1] + last10s['sitk3dl'][:,1]]).T
	last10s['sttk3dr'] = np.array([l, last10s['set3dr'][:,1] + last10s['sitk3dr'][:,1]]).T

	return last10s

#########################################################################################################

def add_ft44_b2(where = ".", last10s = None):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	b2fgmtry = read_b2fgmtry(where = where, verbose = False)

	(fort44, _) = read_ft44(where = where, verbose = False)

	# variables = ['d', 't', 'rflux', 'pflux', 'reflux', 'peflux']
	variables = fort44.keys()
	# species   = ['a', 'm', 'i']
	locations = ['a', 'i', 'l', 'r']
	skip = ['natm', 'nmol', 'nion', 'nfla', 'species', 'atms', 'mols', 'ions']

	b2mn = read_b2mn(where = where)

	for variable in variables:
		if variable not in skip:
			for location in locations:

				data = fort44[variable]

				label = "".join([variable, '3d', location])

				# try/except to avoid error when dealing with
				# fort.44's pden*_int* different format

				try:
					last10 = create_profile(b2mn = b2mn, b2fgmtry = b2fgmtry,
											last10s = last10s, data = data, location = location)


					last10[:,1] = np.abs(last10[:,1])
					last10s[label] = last10

				except: last10s[label] = []

	return last10s

#########################################################################################################

def add_na(where = ".", last10s = None):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	b2fgmtry = read_b2fgmtry(where = where, verbose = False)
	b2fstate = read_b2fstate(where = where, verbose = False)

	locations = ['a', 'i', 'l', 'r']

	b2mn = read_b2mn(where = where)
	data = b2fstate['na']

	for location in locations:

		last10 = create_profile(b2mn = b2mn, b2fgmtry = b2fgmtry,
								last10s = last10s, data = data, location = location)
		last10s['na3d' + location] = last10

	return last10s

#########################################################################################################

def add_Zeff(where = ".", last10s = None):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # June 2022

	b2fstate = read_b2fstate(where = where, verbose = False)
	Za = b2fstate['zamin']

	locations = ['a', 'i', 'l', 'r']

	for location in locations:

		try:

			ne = last10s['ne3d' + location][1:-1,:]
			na = last10s['na3d' + location]

			nae  = np.zeros((na.shape))
			Zeff = np.zeros((ne.shape))
			nae[:,0]  = ne[:,0]
			Zeff[:,0] = ne[:,0]

			for i in range(1,na.shape[1]-1):
				nae[:,i] = na[:,i] / ne[:,1]
				Zeff[:,1] += Za[i]**2 * nae[:,i]

			last10s['nae3d' + location]  = nae
			last10s['Zeff3d' + location] = Zeff

		except: pass

	return last10s

#########################################################################################################

def add_fht(where = ".", last10s = None):

	b2fgmtry = read_b2fgmtry(where = where)
	b2fplasmf = read_b2fplasmf(where = where)

	dummy = np.array([0])

	index = 1

	dr = last10s['ft3dr'][:,0]
	ft = b2fplasmf['fht'][-1,:,index] / b2fgmtry['gs'][-1,:,0] / np.abs(b2fgmtry['bb'][-1,:,0] / b2fgmtry['bb'][-1,:,3])
	ft = np.concatenate((dummy, ft, dummy))
	last10s['fpt3dr'] = np.array([dr, ft]).T

	# ix = 1 because flux from left cell is needed
	
	dl = last10s['ft3dl'][:,0]
	ft = b2fplasmf['fht'][1,:,index] / b2fgmtry['gs'][1,:,0] / np.abs(b2fgmtry['bb'][1,:,0] / b2fgmtry['bb'][1,:,3])
	ft = np.concatenate((dummy, ft, dummy))
	last10s['fpt3dl'] = np.array([dl, ft]).T

	return last10s

#########################################################################################################

def add_grazing_angle(where = ".", last10s = None):

	# for key in ['ft3dlP', 'ft3drP', 'ft3dtlP', 'ft3dtrP']:
	# 	try:
	# 		ga = np.abs(np.arcsin(last10s[key[:-1]][:,1] / last10s[key][:,1]) / np.pi * 180) # [deg]
	# 		last10s['ga' + key[2:5]] = np.array([last10s[key][:,0], ga]).T
	# 	except:
	# 		pass

	b2fgmtry = read_b2fgmtry(where = where)

	# grazing angle = bx / bb

	RAD2DEG = 180 / np.pi

	raise ValueError('To be fixed! Not coherent with ft3d* and ft3d*P...')

	if b2fgmtry['nncut'] > 0:

		dl = last10s['ft3dl'][1:-1,0]
		gal = np.abs(b2fgmtry['bb'][0,:,0] / b2fgmtry['bb'][0,:,3] * RAD2DEG)

		dr = last10s['ft3dr'][1:-1,0]
		gar = np.abs(b2fgmtry['bb'][-1,:,0] / b2fgmtry['bb'][-1,:,3] * RAD2DEG)
		last10s['ga3dl'] = np.array([dl, gal]).T
		last10s['ga3dr'] = np.array([dr, gar]).T

	else:

		# limiter => just half of the radial cells

		dl = last10s['ft3dl'][:-1,0]
		dr = last10s['ft3dr'][:-1,0]

		i = dl.shape[0]

		gal = np.abs(b2fgmtry['bb'][0,-i:,0] / b2fgmtry['bb'][0,-i:,3] * RAD2DEG)
		gar = np.abs(b2fgmtry['bb'][-1,-i:,0] / b2fgmtry['bb'][-1,-i:,3] * RAD2DEG)
		last10s['ga3dl'] = np.array([dl, gal]).T
		last10s['ga3dr'] = np.array([dr, gar]).T

	return last10s

#########################################################################################################

def add_transport(where = ".", last10s = None):

	try:

		data = np.loadtxt(os.path.join(where, 'da1'))

		i = 0

		########

		diff = []
		go_on = True
		while go_on:
			diff.append(data[i,:])
			go_on = data[i+1,0] < 0
			i += 1
		go_on = True
		while go_on:
			diff.append(data[i,:])
			go_on = data[i+1,0] >= 0
			i += 1
		diff = np.array(diff)

		########

		chi_i = []
		go_on = True
		while go_on:
			chi_i.append(data[i,:])
			go_on = data[i+1,0] < 0
			i += 1
		go_on = True
		while go_on:
			chi_i.append(data[i,:])
			go_on = data[i+1,0] >= 0
			i += 1
		chi_i = np.array(chi_i)

		########

		chi_e = []
		go_on = True
		while go_on:
			chi_e.append(data[i,:])
			go_on = data[i+1,0] < 0
			i += 1
		go_on = True
		while go_on:
			try:
				chi_e.append(data[i,:])
				go_on = data[i+1,0] >= 0
				i += 1
			except: go_on = False
		chi_e = np.array(chi_e)

		########

		last10s['diff']  = diff
		last10s['chi_i'] = chi_i
		last10s['chi_e'] = chi_e

	except:
		pass
	
	return last10s

#########################################################################################################

def add_rates(where = ".", last10s = None):

	b2mn = read_b2mn(where = where)
	b2fgmtry = read_b2fgmtry(where = where)

	try:
		rates = pickle.load(open(os.path.join(where, 'rates.pkl'), 'rb'))
	except:
		rates = compute_rates(where = where)

	for key in rates.keys():

		rt = rates[key]
		base = '_'.join([rt['database'], rt['group'], rt['reaction']])
		r = rt['reaction_rate']

		locations = ['a', 'i', 'l', 'r']

		if r != [] and (rt['database'] == 'AMJUEL' or rt['database'] == 'HYDHEL'):

			for location in locations:

				# data = np.zeros((r.shape[1],1))
				label = "".join([base, '_3d', location])

				last10 = create_profile(b2mn = b2mn, b2fgmtry = b2fgmtry,
										last10s = last10s, data = r, location = location)
				last10s[label] = last10

	return last10s