import os
import time
import pickle
import numpy as np

#########################################################################################################

def find_occurrences(target = None, string = None):

	return np.where(np.array([target == character for characer in string]) == True)[0]

#########################################################################################################

def is_uptodate(where = None, what = None):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	time_pkl = os.path.getmtime(os.path.join(where, what + ".pkl"))
	time_original = os.path.getmtime(os.path.join(where, what))

	if time_original < time_pkl: return True
	return False

#########################################################################################################

def load_pickle(where = None, verbose = True, what = None):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	if is_uptodate(where = where, what = what) is False:
		raise ValueError('Pickle file NOT up-to-date :(')

	if verbose is True: print(what + " loaded from .pkl")

	return pickle.load(open(os.path.join(where, what + ".pkl"), "rb"))

#########################################################################################################

def find_sp44(what = None, fort44 = None):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	if   (what.split('b2')[0])[-1] == 'a': sp44 = fort44['neut']['atms']
	elif (what.split('b2')[0])[-1] == 'm': sp44 = fort44['neut']['mols']
	elif (what.split('b2')[0])[-1] == 'i':
		if   (what.split('b2')[0])[-2] == 'a':
			sp44 = fort44['neut']['atms']
		elif (what.split('b2')[0])[-2] == 'm':
			sp44 = fort44['neut']['mols']
		elif (what.split('b2')[0])[-2] == 'n':
			sp44 = ['']
		else:
			sp44 = fort44['neut']['ions']
	elif (what.split('b2')[0])[-1] == 'e':
		if   (what.split('b2')[0])[-2] == 'a':
			sp44 = fort44['neut']['atms']
		elif (what.split('b2')[0])[-2] == 'm':
			sp44 = fort44['neut']['mols']
		else:
			sp44 = ['']
	elif (what.split('b2')[0])[-1] == 'n': sp44 = ["Total neutral"]
	
	elif what[-1] == 'a': sp44 = fort44['neut']['atms']
	elif what[-1] == 'm': sp44 = fort44['neut']['mols']
	elif what[-1] == 'i': sp44 = fort44['neut']['ions']
	elif what[-1] == 'n': sp44 = ["Total neutral"]

	elif what.rfind('neut') > -1: sp44 = fort44['neut']['atms']
	elif what.rfind('mol') > -1:  sp44 = fort44['neut']['mols']
	elif what.rfind('ml') > -1:   sp44 = fort44['neut']['mols']
	elif what.rfind('ion') > -1:  sp44 = fort44['neut']['ions']

	elif what == 'eneutrad': sp44 = fort44['neut']['atms']
	elif what == 'emiss':    sp44 = fort44['neut']['atms'][0] # H only - asuming no H/D/T mixture

	else: sp44 = ['']	

	return sp44

#########################################################################################################

def populate_b2(what = None, b2mn = None, b2fstate = None, b2fplasmf = None, fort44 = None):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	spb2 = b2fstate['species']

	if b2mn['b2mndr_eirene'] == 1: sp44 = find_sp44(what = what, fort44 = fort44)

	try:
		vs = b2fplasmf[what]
		try: 
			if   what[-1] == 'i' or what[-1] == 'a' or what[:2] == 'rq':
				sp = spb2
			elif what[-1] == 'e':
				if what[-2] == 'a':
					sp = spb2
				else:
					sp = ['e-']
		except: sp = ['']
	except:
		if b2mn['b2mndr_eirene'] == 1:
			try:    vs = fort44['neut'][what]
			except: vs = fort44['wld'][what]
			sp = sp44
		else:
			pass

	if not 'sp' in locals(): return ([''], vs)

	return (sp, vs)

#########################################################################################################

def rearrange_quadrangles(crx = None, cry = None):

	# - reorder cr* elements [0,1,2,3] => [0,1,3,2]
	# - add 0th element in last position to plot close polygon
				
	tmp = np.zeros((crx.shape[0], crx.shape[1], crx.shape[2]+1))
	tmp[:,:,0] = crx[:,:,0]
	tmp[:,:,1] = crx[:,:,1]
	tmp[:,:,2] = crx[:,:,3]
	tmp[:,:,3] = crx[:,:,2]
	tmp[:,:,4] = crx[:,:,0]
	crx = tmp

	tmp = np.zeros((cry.shape[0], cry.shape[1], cry.shape[2]+1))
	tmp[:,:,0] = cry[:,:,0]
	tmp[:,:,1] = cry[:,:,1]
	tmp[:,:,2] = cry[:,:,3]
	tmp[:,:,3] = cry[:,:,2]
	tmp[:,:,4] = cry[:,:,0]
	cry = tmp

	return (crx, cry)

#########################################################################################################

def rearrange_triangles(cells = None):

	# - add 0th element in last position to plot close polygon
				
	tmp = np.zeros((cells.shape[0], cells.shape[1]+1))
	tmp[:,0] = cells[:,0]
	tmp[:,1] = cells[:,1]
	tmp[:,2] = cells[:,2]
	tmp[:,3] = cells[:,0]
	cells = tmp

	return cells

#########################################################################################################

def compute_integral(b2fgmtry = None, last10s = None, what = None):

	if   what[-1] == 'l': target = 0
	elif what[-1] == 'r': target = -1
	else: raise ValueError('Integral value can be computed at the targets only!')

	hz = 2 * np.pi * b2fgmtry['crx'][target,:,:].mean(axis = 1)

	# skip additional cell for compatibility (negligible importance)
	dl = np.abs(last10s[what][:-1,0] - last10s[what][1:,0])
	dl = dl[1:]
	values = last10s[what][1:-1,1]

	integral = hz * dl * values

	return integral.sum()