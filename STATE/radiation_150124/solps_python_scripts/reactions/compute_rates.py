import os
import pickle
import numpy as np
import math
import warnings

from solps_python_scripts.read_b2mn import read_b2mn
from solps_python_scripts.read_b2fgmtry import read_b2fgmtry
from solps_python_scripts.read_b2fstate import read_b2fstate
from solps_python_scripts.read_ft44     import read_ft44
from solps_python_scripts.read_ft46     import read_ft46

from solps_python_scripts.reactions.read_reactions import read_all_reactions

############################################################################

def impose_eirene_limits(b2mn = None, n = None, te = None, ti = None):

	# min/max B2 density/temperature used in EIRENE reaction rates
	lims = {}

	for variable in ['na', 'te', 'ti']:
		for extremum in ['min', 'max']:
			try:
				lims['_'.join([variable, extremum])] = float(b2mn['_'.join(['eirene', variable, extremum])])
				if variable == 'na':
					try:
						n[n < lims['na_min']] = lims['na_min']
						print()
						print('    eirene_na_min violated...')
						print()
					except: pass
					try:
						n[n > lims['na_max']] = lims['na_max']
						print()
						print('    eirene_na_max violated...')
						print()
					except: pass
				elif variable == 'te':
					try:
						te[te < lims['te_min']] = lims['te_min']
						print()
						print('    eirene_te_min violated...')
						print()
					except: pass
					try:
						te[te > lims['te_max']] = lims['te_max']
						print()
						print('    eirene_te_max violated...')
						print()
					except: pass
				elif variable == 'ti':
					try: 
						ti[ti < lims['ti_min']] = lims['ti_min']
						print()
						print('    eirene_ti_min violated...')
						print()
					except: pass
					try:
						ti[ti > lims['ti_max']] = lims['ti_max']
						print()
						print('    eirene_ti_max violated...')
						print()
					except: pass
			except: pass

	return (n, te, ti)

############################################################################

def lnsv_fit(coef = None, v0 = None, v1 = None):

	lnsv = np.zeros(v0.shape, dtype = np.float128)

	if v1 is None:
		here = v0 > 0
		for n in range(coef.shape[0]):
			lnsv[here] += coef[n] * (np.log(v0[here]))**n

	elif v1 is not None:
		here = (v0 > 0) * (v1 > 0)
		for n in range(coef.shape[0]):
			for m in range(coef.shape[1]):
				try:
					lnsv[here] += coef[n,m] * (np.log(v0[here]))**n * (np.log(v1[here]))**m
				except:
					lnsv[here] += coef[n,m] * (np.log(v0[here]))**n * (np.log(v1))**m

	return lnsv

############################################################################

def evaluate_rate(reaction = None, fort44 = None, n0 = None, t0 = None, vth0 = None, n1 = None, t1 = None, vth1 = None):

	# unit conversion needed:
	# - 1E-06 from [cm^3] to [m^3]
	# - 1E-08 AMJUEL/HYDHEL convention
	factor = 1E-06 * 1E-08

	# lnsv_plot to plot fit (not reaction rate)
	#
	# Convention isinstance(lnsv_plot, dict):
	# - if lnsv_plot.keys() == 'all' => one-parameter fit
	# - if lnsv_plot.keys() != 'all' => two-parameter fit (2nd parameter = key)

	num0 = int(1E+02)
	num1 = int(5E+00)
	minimum = np.log10(t0[t0 > 0].min())
	maximum = np.log10(t0.max())
	v0_plot = np.logspace(minimum, maximum, num0)
	v1_plot = np.ones(num1) * np.nan

	database = reaction['database']
	group    = reaction['group']
	coef     = reaction['fit_coefficients']

	if database == 'AMJUEL' or database == 'HYDHEL': 

		if group == 'H.2':
			lnsv = lnsv_fit(coef = coef, v0 = t0, v1 = None) 
			lnsv_plot = lnsv_fit(coef = coef, v0 = v0_plot, v1 = None)
			lnsv_plot = (v0_plot, v1_plot, lnsv_plot)

		elif group == 'H.8':
			# see AMJUEL manual pag. 247: from lnsvE to lnsv
			lnsvE = lnsv_fit(coef = coef, v0 = t0, v1 = None) 
			lnsv = np.log(np.exp(lnsvE) / 8.96E-01 / t0)	
			lnsv_plot = lnsv_fit(coef = coef, v0 = v0_plot, v1 = None)
			lnsv_plot = (v0_plot, v1_plot, np.log(np.exp(lnsv_plot) / 8.96E-01 / v0_plot))
			
		elif group == 'H.3':
			if 'H' not in fort44['atms'] and 'D' not in fort44['atms']:
				raise ValueError('H/D must be present in H.3-type reaction!')
			else:
				if 'H_2' in reaction['species_0'] or 'H_2' in reaction['species_1']:
					e0 = fort44['e0mb2']
				elif 'H' in reaction['species_0'] or 'H' in reaction['species_1']:
					e0 = fort44['e0ab2']
			# lnsv = lnsv_fit(coef = coef, v0 = t0, v1 = t1) 
			lnsv = lnsv_fit(coef = coef, v0 = t0, v1 = e0) 
			# 2nd tuple element => e0
			dummy_plot = np.zeros((num0, num1))
			e0_plot = np.logspace(np.log10(e0[e0 > 0].min()), np.log10(e0.max()), num1)
			for n in range(num1):
				dummy_plot[:,n] = lnsv_fit(coef = coef, v0 = v0_plot, v1 = e0_plot[n])
			lnsv_plot = (v0_plot, e0_plot, dummy_plot)

		elif group == 'H.4' or group == 'H.12':
			lnsv = lnsv_fit(coef = coef, v0 = t0, v1 = n0 * factor) 
			# 2nd tuple element => n0
			dummy_plot = np.zeros((num0, num1))
			n0_plot = np.logspace(np.log10(n0[n0 > 0].min()), np.log10(n0.max()), num1)
			for n in range(num1):
				dummy_plot[:,n] = lnsv_fit(coef = coef, v0 = v0_plot, v1 = n0_plot[n] * factor)
			lnsv_plot = (v0_plot, n0_plot, dummy_plot)

		else:
			print()
			print('  H.0, H.1, H.10 not yet implemented!')
			# lnsv = -inf is returned => exp(lnsv) = 0
			lnsv = np.ones(n0.shape) * (-1) * math.inf
			lnsv_plot = []

	elif database == 'AMMONX':
		# CAUTION.
		#
		# - neutral gas temperature needed
		# - temperature assumed is sp1's (sp0 moving in background sp1)
		#
		# => R-H-H2 != R-H2-H because tdena != tdenm
		v0 = t1
		lnsv = lnsv_fit(coef = coef, v0 = v0, v1 = None)
		minimum = np.log10(t1[t1 > 0].min())
		maximum = np.log10(t1.max())
		v0_plot = np.logspace(minimum, maximum, num0)
		lnsv_plot = lnsv_fit(coef = coef, v0 = v0_plot, v1 = None)
		lnsv_plot = (v0_plot, v1_plot, lnsv_plot)

	else:
		raise ValueError(database + ' not yet implemented!')

	# reaction rate

	if   group != 'H.12': R01 = n0 * n1 * np.exp(lnsv) * 1E-06
	elif group == 'H.12':
		if   reaction['reaction'] == '2.1.5tot': R01 = lnsv # different convention
		elif reaction['reaction'] != '2.1.5tot': R01 = np.exp(lnsv)

	# mean free path

	here = (R01 > 0) * (np.isnan(R01) == 0)
	l0 = np.zeros(lnsv.shape)
	l1 = np.zeros(lnsv.shape)

	if database == 'AMJUEL' or database == 'HYDHEL': 
		
		if vth0 is not None: l0[here] = vth0[here] * n0[here] / R01[here]
		if vth1 is not None: l1[here] = vth1[here] * n1[here] / R01[here]

	elif database == 'AMMONX':

		sigma = 1.5E-19 # [m^{-2}] from EIRENE manual
		l0 = 1 / (np.sqrt(2) * sigma * n1)
		l1 = 1 / (np.sqrt(2) * sigma * n0)

	l0[np.isnan(l0)] = 0
	l1[np.isnan(l1)] = 0
	l0[np.isinf(l0)] = 0
	l1[np.isinf(l1)] = 0

	# 1E-06 to go from [s^-1 cm^-3] to [s^-1 m^-3]
	return (R01, l0, l1, lnsv_plot)

############################################################################

def compute_rates(where = ".", auto = True, verbose = True, save = True):

	# CAUTION.
	#
	# - units: [cm^3 / s]
    # 
    # - H.1 => potential => NEGLECTED
    # - H.2 => ion temperature dependence
    # - H.3 => ion & neutral temperature dependence
    # - H.4 => ion density & temperature dependence
    # - H.8 => electron temperature dependence + manual pag. 247 (sigma * v * E)
    # - H.10 => energy loss => NEGLECTED

	if auto is True:

		try:

			rates = pickle.load(open(os.path.join(where, 'rates.pkl'), 'rb'))
			if verbose is True: print('compute_rates - rates.pkl')
			return rates

		except:

			if verbose is True: print('compute_rates ...')

			b2mn        = read_b2mn(where = where)
			b2fgmtry    = read_b2fgmtry(where = where, verbose = True)
			b2fstate    = read_b2fstate(where = where, verbose = True)
			(fort44, _) = read_ft44(where = where, verbose = True)
			fort46      = read_ft46(where = where, verbose = True)

			reactions = read_all_reactions(where = where)
			rates = reactions

			for key in reactions.keys():

				# CAUTION.
				#
				# - species index following B2 convention for plasma species
				# - electrons, if any, always first 
				# - for AMMONX ~ BGK then 1 is neutral atoms

				reaction = reactions[key]
				database = reaction['database']

				if (database == 'AMJUEL' or database == 'HYDHEL') and reaction['fit_coefficients'] != []:

					n = np.zeros((b2fgmtry['nx'], b2fgmtry['ny'], 2))
					t = np.zeros((b2fgmtry['nx'], b2fgmtry['ny'], 2))
					vth = np.zeros((b2fgmtry['nx'], b2fgmtry['ny'], 2))

					sp0 = reaction['species_0']
					sp1 = reaction['species_1']

					for sp in [sp0, sp1]:

						# homogenise fort.44 species notation with database LaTeX's

						isp = 1

						# electrons: always isp = 0
						if sp == 'e':
							isp = 0
							n[:,:,isp] = b2fstate['ne']
							t[:,:,isp] = b2fstate['te']
							vth[:,:,isp] = np.sqrt(b2fstate['te'] * 1.6E-19 / 9E-31)
							(n[:,:,isp], t[:,:,isp], _) = impose_eirene_limits(b2mn = b2mn,
																	  n  = n[:,:,isp],
																	  te = t[:,:,isp],
																	  ti = None)
						# main ions: always isp = 0, unless electrons are present => isp = 1
						elif sp == 'p' or sp == 'H^+':
							if 'e' not in [sp0, sp1]: isp = 0
							n[:,:,isp] = b2fstate['na'][:,:,1]
							t[:,:,isp] = b2fstate['ti']
							vth[:,:,isp] = np.sqrt(b2fstate['ti'] * 1.6E-19 / 1.67E-27 / b2fstate['am'][1])
							(n[:,:,isp], _, t[:,:,isp]) = impose_eirene_limits(b2mn = b2mn,
																	  n  = n[:,:,isp],
																	  te = None,
																	  ti = t[:,:,isp])
						# atoms
						elif sp in fort44['atms']:
							n[:,:,isp] = fort44['dab2'][:,:, fort44['atms'].index(sp)]
							t[:,:,isp] = fort44['tab2'][:,:, fort44['atms'].index(sp)]
							vth[:,:,isp] = fort44['vtab2'][:,:, fort44['atms'].index(sp)]
						# molecules
						elif sp in fort44['mols']:
							n[:,:,isp] = fort44['dmb2'][:,:, fort44['mols'].index(sp)]
							t[:,:,isp] = fort44['tmb2'][:,:, fort44['mols'].index(sp)]
							vth[:,:,isp] = fort44['vtmb2'][:,:, fort44['mols'].index(sp)]
						# molecular ions
						elif sp in fort44['ions']:
							n[:,:,isp] = fort44['dib2'][:,:, fort44['ions'].index(sp)]
							t[:,:,isp] = fort44['tib2'][:,:, fort44['ions'].index(sp)]	
							# vth[:,:,isp] just zeros
						# main neutral atoms
						elif sp == 'H' or sp == 'H(1s)':
							n[:,:,isp] = fort44['dab2'][:,:,0]
							t[:,:,isp] = fort44['tab2'][:,:,0]
							vth[:,:,isp] = fort44['vtab2'][:,:,0]
						# main neutral molecules
						elif sp == 'H_2' or sp == 'D_2':
							n[:,:,isp] = fort44['dmb2'][:,:,0]
							t[:,:,isp] = fort44['tmb2'][:,:,0]
							vth[:,:,isp] = fort44['vtmb2'][:,:,0]
						# main neutral molecular ions
						elif sp == 'H_2^+' or sp == 'D_2^+' or sp == 'H_2^+(v)':
							n[:,:,isp] = fort44['dib2'][:,:,0]
							t[:,:,isp] = fort44['tib2'][:,:,0]
							# vth[:,:,isp] just zeros
						else:
							print()
							print('  Species {} switch not implemented!'.format(sp))

					# removing nans
					n[np.isnan(n)] = 0
					t[np.isnan(t)] = 0

					# raise ValueError('Fix which density/temperature is used in the fit!')

					(rate, mfp0, mfp1, reactions[key]['lnsv']) = evaluate_rate(reaction = reaction, fort44 = fort44,
															                   n0 = n[:,:,0], t0 = t[:,:,0], vth0 = vth[:,:,0],
															                   n1 = n[:,:,1], t1 = t[:,:,1], vth1 = vth[:,:,1])

				elif database == 'AMMONX' and reaction['fit_coefficients'] != []:

					n = np.zeros((fort46['pdena'].shape[0], 2))
					t = np.zeros((fort46['pdena'].shape[0], 2))
					vth = np.zeros((fort46['pdena'].shape[0], 2))

					# convention: R-sp0-sp1, e.g. R-H-H2

					sp0 = reaction['reaction'].split('-')[1]
					sp1 = reaction['reaction'].split('-')[2]

					for (isp, sp) in enumerate([sp0, sp1]):

						# homogenise fort.44 species notation with database LaTeX's

						# main neutral atoms
						if sp == 'H':
							n[:,isp] = fort46['pdena'][:,0]
							t[:,isp] = fort46['tdena'][:,0]
						# main neutral molecules
						elif sp == 'H2':
							n[:,isp] = fort46['pdenm'][:,0]
							t[:,isp] = fort46['tdenm'][:,0]
						else:
							print()
							print('  Species {} switch not implemented!'.format(sp))

					# removing nans
					n[np.isnan(n)] = 0
					t[np.isnan(t)] = 0

					(rate, mfp0, mfp1, reactions[key]['lnsv']) = evaluate_rate(reaction = reaction,
																		   	   n0 = n[:,0], t0 = t[:,0], vth0 = vth[:,0],
																		   	   n1 = n[:,1], t1 = t[:,1], vth1 = vth[:,1])

				else:
					rate = []
					mfp0 = []
					mfp1 = []
					reactions[key]['lnsv'] = []

				rates[key]['reaction_rate'] = rate
				rates[key]['mfp_0'] = mfp0
				rates[key]['mfp_1'] = mfp1

			if save is True:
				pickle.dump(rates, open(os.path.join(where, 'rates.pkl'), 'wb'))

	else:
		raise ValueError('auto = False not yet implemented!')

	return rates