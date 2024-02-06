import os
import csv
import numpy as np
import matplotlib.pyplot as plt

from decimal import Decimal

from raysect.core import Point2D, Point3D, translate, rotate_z
from raysect.optical.material import VolumeTransform
from raysect.primitive import Cylinder, Subtract

from cherab.tools.equilibrium import import_eqdsk
from cherab.core.math import Interpolate1DLinear, Interpolate2DLinear
from cherab.core.math.mappers import AxisymmetricMapper, DiscreteToroidalMapper
from cherab.tools.emitters import RadiationFunction

from quantum_data import FractionExcitedHydrogenGas

# from create_step_function import CreateStepFunction
# from create_absorption_function import CreateAbsorptionFunction

from solps_python_scripts.reactions.read_reactions import read_one_reaction
from solps_python_scripts.reactions.compute_rates  import evaluate_rate

from accessories import plot_1d, plot_2d, ComputeDecayLength, make_chamber_wall
from accessories import normalise_psi, ExponentialDecay, Gaussian, FermiDirac, DensityDoublingStangeby

##############################################################################################
##############################################################################################
##############################################################################################

qe = 1.6E-19
me  = 9.0E-31
mD  = 2 * 1.67E-27
mD2 = 2 * mD

##############################################################################################

def make_homemade_distributions(cfg = None, verbose = False):

	location = os.path.join(cfg['baserun'],
							'input',
							cfg['run'],
							cfg['plasma']['ASTRA']['ASTRA_directory'])

	eqFile = cfg['plasma']['ASTRA']['geqdsk_name']
	eqTime = cfg['plasma']['ASTRA']['geqdsk_time']

	equilibrium = import_eqdsk(os.path.join(location, eqFile + '_' + eqTime + '.geqdsk'))

	# equilibrium = import_eqdsk('/home/matteo.moscheni/raysect-cherab/runs/h_alpha_camera/input/9229_equ@65ms_HW/ASTRA_data/ST40_9229_EFIT_BEST_65p89ms.geqdsk')

	psiN = normalise_psi(equilibrium)
	psiN_lcfs = 1.0
	psiNmax = psiN.max()
	psiNmin = psiN.min()

	psiN_1d = np.linspace(psiNmin, psiNmax, 1000)

	core = psiN_1d <= 1.0
	sol  = psiN_1d >= 1.0

	ng_parameters      = cfg['plasma']['homemade']['neutral_density']
	ng_sol             = ng_parameters['ng_sol'] # [m^{-3}]
	ng_core            = ng_parameters['ng_core'] # [m^{-3}]
	smoothing_ng_FD    = ng_parameters['smoothing_ng_FD'] # [-]
	smoothing_ng_G     = ng_parameters['smoothing_ng_G'] # [-]
	radial_reduction_G = ng_parameters['radial_reduction_G'] # [-]
	peak_G             = ng_parameters['peak_G'] # [m^{-3}]
	r0_G               = ng_parameters['r0_G'] # [m]
	z01_G              = ng_parameters['z01_G'] # [m]
	z02_G              = ng_parameters['z02_G'] # [m]
	ng_raw = (ng_sol - ng_core) * FermiDirac(psi = psiN_1d, psi_lcfs = psiN_lcfs, smoothing = smoothing_ng_FD) + ng_core
	ng_f1d = Interpolate1DLinear(psiN_1d, ng_raw, extrapolate = True)

	Tp_parameters = cfg['plasma']['homemade']['plasma_temperature']
	Tpmax         = Tp_parameters['Tpmax'] # [eV]
	Tpsep         = Tp_parameters['Tpsep'] # [eV]
	Tpmin         = Tp_parameters['Tpmin'] # [eV]
	ltpsi         = Tp_parameters['ltpsi'] # [-] 
	ltrho         = Tp_parameters['ltrho'] # [m] 
	smoothing_Tp = 1 / np.log(Tpmax / Tpsep)
	Tp_raw_core = Tpmax * Gaussian(x = psiN_1d, x0 = 0.0, smoothing = smoothing_Tp) * core
	Tp_raw_sol  = Tpsep * ExponentialDecay(x = psiN_1d, x0 = 1.0, decay_length = ltpsi) * sol
	Tp_raw = Tp_raw_core + Tp_raw_sol
	Tp_raw[np.where(Tp_raw < Tpmin)] = Tpmin # minimum temperature due to neutral-plasma thermal equilibration
	Tp_f1d = Interpolate1DLinear(psiN_1d, Tp_raw, extrapolate = True)

	# - for Tg and ng definition see [M. Moscheni et al 2022 Nucl. Fusion 62 056009]
	# - Tg = Ti in SOL due to neutral-plasma thermal equilibration 
	# - Tg = Tsep in core due to large mfp (no temperature-equilibrating collisions in small ST40 ~ only ionisations in core)
	Tg_raw = Tpsep * core + Tp_raw * sol
	Tg_f1d = Interpolate1DLinear(psiN_1d, Tg_raw, extrapolate = True)

	np_parameters = cfg['plasma']['homemade']['plasma_density']
	npmax         = np_parameters['npmax'] # [m^{-3}]
	npsep         = np_parameters['npsep'] # [m^{-3}]
	npmin         = np_parameters['npmin'] # [m^{-3}] - needed to use AMJUEL data
	lnpsi         = np_parameters['lnpsi'] # [-] 
	lnrho         = np_parameters['lnrho'] # [m] 
	smoothing_np = 1 / np.log(npmax / npsep)
	np_raw_core = npmax * Gaussian(x = psiN_1d, x0 = 0.0, smoothing = smoothing_np) * core
	np_raw_sol  = npsep * ExponentialDecay(x = psiN_1d, x0 = 1.0, decay_length = lnpsi) * sol
	np_raw = np_raw_core + np_raw_sol
	np_raw[np.where(np_raw < npmin)] = npmin
	np_f1d = Interpolate1DLinear(psiN_1d, np_raw, extrapolate = True)
	
	r = equilibrium.r_data
	z = equilibrium.z_data

	ng_raw_2d = np.zeros((np.size(z), np.size(r)))
	Tg_raw_2d = np.zeros((np.size(z), np.size(r)))
	Tp_raw_2d = np.zeros((np.size(z), np.size(r)))
	np_raw_2d = np.zeros((np.size(z), np.size(r)))
	for i in range(ng_raw_2d.shape[0]):
		for j in range(ng_raw_2d.shape[1]):

			point = Point2D(r[j], z[i])
			rho = equilibrium.magnetic_axis.distance_to(point)
			
			ng_raw_2d[i,j] = ng_f1d(psiN[j,i]) + \
			peak_G * \
			 Gaussian(x = r[j], x0 = r0_G, smoothing = smoothing_ng_G * radial_reduction_G) * \
			(Gaussian(x = z[i], x0 = z01_G, smoothing = smoothing_ng_G) + Gaussian(x = z[i], x0 = z02_G, smoothing = smoothing_ng_G))

			Tg_raw_2d[i,j] = Tg_f1d(psiN[j,i])
			if equilibrium.inside_lcfs(r[j],z[i]) == 0: Tg_raw_2d[i,j] *= ExponentialDecay(x = rho, x0 = 0.0, decay_length = ltrho)
			Tp_raw_2d[i,j] = Tp_f1d(psiN[j,i])
			if equilibrium.inside_lcfs(r[j],z[i]) == 0: Tp_raw_2d[i,j] *= ExponentialDecay(x = rho, x0 = 0.0, decay_length = ltrho)
			np_raw_2d[i,j] = np_f1d(psiN[j,i])
			if cfg['plasma']['homemade']['doubleDensity_2ptModel'] is True:
				np_raw_2d[i,j] *= DensityDoublingStangeby(x = r[j], y = z[i], eq = equilibrium)
			if equilibrium.inside_lcfs(r[j],z[i]) == 0: np_raw_2d[i,j] *= ExponentialDecay(x = rho, x0 = 0.0, decay_length = lnrho)

	Tp_raw_2d[np.where(Tp_raw_2d < Tpmin)] = Tpmin
	np_raw_2d[np.where(np_raw_2d < npmin)] = npmin

	ng_f2d = Interpolate2DLinear(r, z, ng_raw_2d.T, extrapolate = True)
	Tg_f2d = Interpolate2DLinear(r, z, Tg_raw_2d.T, extrapolate = True)
	Tp_f2d = Interpolate2DLinear(r, z, Tp_raw_2d.T,  extrapolate = True)
	np_f2d = Interpolate2DLinear(r, z, np_raw_2d.T, extrapolate = True)

	if cfg["plotting"]["plot_homemade_emission"] is True and verbose is True:

		plot_2d(ri = r, zi = z, f2d = ng_f2d, scale = "log", eq = equilibrium, title = r'Neutral density $[m^{-3}]$', cmap_label = '$log_{10}$')
		plot_2d(ri = r, zi = z, f2d = Tg_f2d, scale = "log", eq = equilibrium, title = r'Neutral temperature $[eV]$')
		plot_2d(ri = r, zi = z, f2d = np_f2d, scale = "log", eq = equilibrium, title = r'Plasma density $[m^{-3}]$', cmap_label = '$log_{10}$')
		plot_2d(ri = r, zi = z, f2d = Tp_f2d, scale = "log", eq = equilibrium, title = r'Plasma temperature $[eV]$')

		plot_1d(x = r, y = ng_f2d, y0 = equilibrium.magnetic_axis.y, scale = "lin",
			    xlabel = '$R \\; [m]$',
			    ylabel = '$n_g = n_{D} + n_{D_2} \\; [m^{-3}]$',
			    title = 'Neutral density (mid-height)')
		plot_1d(x = r, y = Tg_f2d, y0 = equilibrium.magnetic_axis.y, scale = "lin",
			    xlabel = '$R \\; [m]$',
			    ylabel = '$T_g \\; [eV]$',
			    title = 'Neutral temperature (mid-height)')
		(ri, Tp_) = plot_1d(x = r, y = Tp_f2d, y0 = equilibrium.magnetic_axis.y, scale = "lin",
			    xlabel = '$R \\; [m]$',
			    ylabel = '$T_i = T_e \\; [eV]$',
			    title = 'Plasma temperature (mid-height)',
			    output = True)
		(ri, np_) = plot_1d(x = r, y = np_f2d, y0 = equilibrium.magnetic_axis.y, scale = "lin",
			    xlabel = '$R \\; [m]$',
			    ylabel = '$n_e = n_i \\; [m^{-3}]$',
			    title = 'Plasma density (mid-height)',
			    output = True)

		ComputeDecayLength(R = ri, data = np_, datamin = npmin, title = 'Plasma density:')
		ComputeDecayLength(R = ri, data = Tp_, datamin = Tpmin, title = 'Plasma temperature:')

		plot_1d(x = psiN_1d, y = ng_raw,
			    xlabel = '$\\psi \\; [-]$',
			    ylabel = '$n_g = n_{D} + n_{D_2} \\; [m^{-3}]$',
			    title = 'Neutral density')
		plot_1d(x = psiN_1d, y = Tg_raw,
			    xlabel = '$\\psi \\; [-]$',
			    ylabel = '$T_g \\; [eV]$',
			    title = 'Neutral temperature')
		plot_1d(x = psiN_1d, y = Tp_raw * 1E-03,
			    xlabel = '$\\psi \\; [-]$',
			    ylabel = '$T_i = T_e \\; [keV]$',
			    title = 'Plasma temperature')
		plot_1d(x = psiN_1d, y = np_raw,
			    xlabel = '$\\psi \\; [-]$',
			    ylabel = '$n_e = n_i \\; [m^{-3}]$',
			    title = 'Plasma density')

		plot_1d(x = psiN_1d, y = (np_raw, Tp_raw * 1E-03),
			    xlabel = '$\\psi \\; [-]$',
			    ylabel = ('$n_e = n_i \\; [m^{-3}]$', '$T_e = T_i \\; [keV]$'),
			    title = 'Plasma density and temperature')

	return ((ng_raw_2d, ng_f2d), (Tg_raw_2d, Tg_f2d), (np_raw_2d, np_f2d), (Tp_raw_2d, Tp_f2d))

####################################################################################################################################################

def make_homemade_emission(cfg = None, verbose = False):

	location = os.path.join(cfg['baserun'],
							'input',
							cfg['run'],
							cfg['plasma']['ASTRA']['ASTRA_directory'])

	eqFile = cfg['plasma']['ASTRA']['geqdsk_name']
	eqTime = cfg['plasma']['ASTRA']['geqdsk_time']

	equilibrium = import_eqdsk(os.path.join(location, eqFile + '_' + eqTime + '.geqdsk'))

	# equilibrium = import_eqdsk('/home/matteo.moscheni/raysect-cherab/runs/h_alpha_camera/input/9229_equ@65ms_HW/ASTRA_data/ST40_9229_EFIT_BEST_65p89ms.geqdsk')

	(ng_, Tg_, np_, Tp_) = make_homemade_distributions(cfg = cfg, verbose = verbose)

	(ng_raw_2d, ng_f2d) = ng_
	(Tg_raw_2d, Tg_f2d) = Tg_
	(np_raw_2d, np_f2d) = np_
	(Tp_raw_2d, Tp_f2d) = Tp_

	# Boltzmann-like (i.e. neutral gas) fraction of excited hydrogen
	fraction_3_Boltzmann = FractionExcitedHydrogenGas(n_quantum = 3, T_eV = Tg_raw_2d)
	# sim.Ha_emission = sim.Ha_Einstein_coeff * f_3 * sim.pdena

	ri = equilibrium.r_data
	zi = equilibrium.z_data

	coef_21   = read_one_reaction(database = 'AMJUEL', group = 'H.12', reaction = '2.1.5b')
	coef_31   = read_one_reaction(database = 'AMJUEL', group = 'H.12', reaction = '2.1.5a')
	coef_tot1 = read_one_reaction(database = 'AMJUEL', group = 'H.12', reaction = '2.1.5tot')

	reaction_21 = {
		'database': 'AMJUEL',
		'group':    'H.12',
		'reaction': '2.1.5b',
		'fit_coefficients': coef_21
	}
	reaction_31 = {
		'database': 'AMJUEL',
		'group':    'H.12',
		'reaction': '2.1.5a',
		'fit_coefficients': coef_31
	}
	reaction_tot1 = {
		'database': 'AMJUEL',
		'group':    'H.12',
		'reaction': '2.1.5tot',
		'fit_coefficients': coef_tot1
	}

	(fraction_21_Amjuel, _, _, lnsv_21)     = evaluate_rate(reaction = reaction_21,   n0 = np_raw_2d, t0 = Tp_raw_2d, n1 = ng_raw_2d, t1 = Tg_raw_2d)
	(fraction_31_Amjuel, _, _, lnsv_31)     = evaluate_rate(reaction = reaction_31,   n0 = np_raw_2d, t0 = Tp_raw_2d, n1 = ng_raw_2d, t1 = Tg_raw_2d)
	(fraction_tot1_Amjuel, _, _, lnsv_tot1) = evaluate_rate(reaction = reaction_tot1, n0 = np_raw_2d, t0 = Tp_raw_2d, n1 = ng_raw_2d, t1 = Tg_raw_2d)

	fraction_2_Amjuel = fraction_21_Amjuel / fraction_tot1_Amjuel
	fraction_3_Amjuel = fraction_31_Amjuel / fraction_tot1_Amjuel

	# [s^{-1}]: Einstein coefficient for H*(3)->H*(2) transition
	Da_Einstein_coeff = 4.41E+07

	if cfg['plasma']['homemade']['sum_Amjuel_and_Boltzmann'] is True:
		fraction_2_total = fraction_2_Amjuel + fraction_2_Boltzmann
		fraction_3_total = fraction_3_Amjuel + fraction_3_Boltzmann
	else:
		fraction_2_total = fraction_2_Amjuel
		fraction_3_total = fraction_3_Amjuel

	ng_2_raw_2d = ng_raw_2d * fraction_2_total
	ng_3_raw_2d = ng_raw_2d * fraction_3_total

	Da_raw_2d = ng_3_raw_2d * Da_Einstein_coeff
	# Da_raw_2d[np.where(Da_raw_2d < 1e18)] = 1e18
	# Da_raw_2d[np.where(Da_raw_2d > 6e21)] = 6e21

	ng_2_f2d = Interpolate2DLinear(ri, zi, ng_2_raw_2d.T, extrapolate = True)
	ng_3_f2d = Interpolate2DLinear(ri, zi, ng_3_raw_2d.T, extrapolate = True)
	Da_f2d = Interpolate2DLinear(ri, zi, Da_raw_2d.T, extrapolate = True)
	f3_f2d = Interpolate2DLinear(ri, zi, fraction_3_total.T, extrapolate = True)
	f2_f2d = Interpolate2DLinear(ri, zi, fraction_2_total.T, extrapolate = True)

	if cfg['plotting']['plot_homemade_emission'] is True and verbose is True:

		i = 0
		titles = ['$D^*(3)/D^*(1)$', '$D_{tot}/D^*(1)$']
		left   = [1e-1,  1e0]
		right  = [1e3,   2e4]
		top    = [1e-1,  1.1]
		bottom = [1e-12, 1.0]
		ylabel  = ['$[-]$', '$[-]$']

		for lnsv in (lnsv_31, lnsv_tot1):
			(p0, p1, sv) = lnsv
			if i != 1: sv = np.exp(sv)
			# elif i == 1: sv = 10**(sv) # different IMPLICIT DAMN convention: sv is already what is needed 
			fig, ax = plt.subplots()
			for j in range(p1.shape[0]):
				parameter = '$n =$ '
				units     = ' $cm^{-3}$' 
				if i != 1: ax.loglog(p0, sv[:,j],   label = parameter + '%.1E' % Decimal(str(p1[j] * 1E-06)) + units)
				else:      ax.semilogx(p0, sv[:,j], label = parameter + '%.1E' % Decimal(str(p1[j] * 1E-06)) + units)
			plt.xlabel(r'$T_i = T_e \; [eV]$')
			plt.ylabel(ylabel[i])
			plt.xlim(left = left[i], right = right[i])
			plt.ylim(top = top[i], bottom = bottom[i])
			plt.title(titles[i])
			ax.legend(loc = 'best')
			plt.grid(True)
			i += 1

		if cfg['plasma']['homemade']['sum_Amjuel_and_Boltzmann'] is True:
			plot_2d(ri = ri, zi = zi, f2d = fraction_3_Boltzmann, eq = equilibrium, title = r'Boltzmann-like $D^*(3)/D_{tot}$ fraction')
			plot_2d(ri = ri, zi = zi, f2d = fraction_3_Amjuel, eq = equilibrium, title = r'Amjuel-like $D^*(3)/D_{tot}$ fraction')
		plot_2d(ri = ri, zi = zi, f2d = f3_f2d * 1E+02, eq = equilibrium, title = r'$D^*(3)$ fractional abundance [%]')
		plot_2d(ri = ri, zi = zi, f2d = f2_f2d * 1E+02, eq = equilibrium, title = r'$D^*(2)$ fractional abundance [%]')

		_ = compute_mfp_ionz(cfg = cfg, verbose = True)

	return ((ng_2_f2d, ng_2_raw_2d), (ng_3_f2d, ng_3_raw_2d), (Da_f2d, Da_raw_2d))

####################################################################################################################################################

def compute_mfp_ionz(cfg = None, verbose = False):

	location = os.path.join(cfg['baserun'],
							'input',
							cfg['run'],
							cfg['plasma']['ASTRA']['ASTRA_directory'])

	eqFile = cfg['plasma']['ASTRA']['geqdsk_name']
	eqTime = cfg['plasma']['ASTRA']['geqdsk_time']

	equilibrium = import_eqdsk(os.path.join(location, eqFile + '_' + eqTime + '.geqdsk'))

	(ng_, Tg_, np_, Tp_) = make_homemade_distributions(cfg = cfg, verbose = False)

	(ng_raw_2d, ng_f2d) = ng_
	(Tg_raw_2d, Tg_f2d) = Tg_
	(np_raw_2d, np_f2d) = np_
	(Tp_raw_2d, Tp_f2d) = Tp_

	ri = equilibrium.r_data
	zi = equilibrium.z_data

	coef_ionz_atm   = read_one_reaction(database = 'AMJUEL', group = 'H.4', reaction = '2.1.5')
	coef_ionz_mol   = read_one_reaction(database = 'AMJUEL', group = 'H.4', reaction = '2.2.9')

	reaction_ionz_atm = {
		'database': 'AMJUEL',
		'group':    'H.4',
		'reaction': '2.1.5',
		'fit_coefficients': coef_ionz_atm
	}
	reaction_ionz_mol = {
		'database': 'AMJUEL',
		'group':    'H.4',
		'reaction': '2.2.9',
		'fit_coefficients': coef_ionz_mol
	}

	(rate_ionz_atm, _, mfp1_ionz_atm, lnsv_ionz_atm) = \
		evaluate_rate(reaction = reaction_ionz_atm,
					  n0 = np_raw_2d, t0 = Tp_raw_2d, vth0 = np.sqrt(2 * qe * Tp_raw_2d / me),
					  n1 = ng_raw_2d, t1 = Tg_raw_2d, vth1 = np.sqrt(2 * qe * Tg_raw_2d / mD))
	mfp1_ionz_atm[np.where(mfp1_ionz_atm > 1E+03)] = 1E+03

	(rate_ionz_mol, _, mfp1_ionz_mol, lnsv_ionz_mol) = \
		evaluate_rate(reaction = reaction_ionz_mol,
					  n0 = np_raw_2d, t0 = Tp_raw_2d, vth0 = np.sqrt(2 * qe * Tp_raw_2d / me),
					  n1 = ng_raw_2d, t1 = Tg_raw_2d, vth1 = np.sqrt(2 * qe * Tg_raw_2d / mD2))
	mfp1_ionz_mol[np.where(mfp1_ionz_mol > 1E+03)] = 1E+03

	total_rate_f2d = rate_ionz_atm + rate_ionz_mol

	mfp = np.minimum(mfp1_ionz_atm, mfp1_ionz_mol)
	mfp_f2d = Interpolate2DLinear(ri, zi, mfp.T, extrapolate = True)
	total_rate_f2d = Interpolate2DLinear(ri, zi, total_rate_f2d.T, extrapolate = True)

	if cfg['plotting']['plot_homemade_emission'] is True and verbose is True:

		i = 0
		titles = ['Atomic ionization', 'Molecular ionization']
		left   = [1e-1, 1e-1]
		right  = [1e4, 1e3]
		top    = [1e-6, 1e-7]
		bottom = [1e-14, 1e-12]

		# for lnsv in (lnsv_ionz_atm, lnsv_ionz_mol):
		# 	(p0, p1, sv) = lnsv
		# 	sv = np.exp(sv) # [cm^3/s]
		# 	fig, ax = plt.subplots()
		# 	for j in range(p1.shape[0]):
		# 		parameter = '$n =$ '
		# 		units     = ' $cm^{-3}$' 
		# 		ax.loglog(p0, sv[:,j],   label = parameter + '%.1E' % Decimal(str(p1[j] * 1E-06)) + units)
		# 	plt.xlabel(r'$T_i = T_e \; [eV]$')
		# 	plt.ylabel('$[cm^3 \\cdot s^{-1}]$')
		# 	plt.xlim(left = left[i], right = right[i])
		# 	plt.ylim(top = top[i], bottom = bottom[i])
		# 	plt.title(titles[i])
		# 	ax.legend(loc = 'best')
		# 	plt.grid(True)
		# 	i += 1

		plot_1d(x = ri, y = mfp_f2d, y0 = equilibrium.magnetic_axis.y, scale = "log",
			    xlabel = '$R \\; [m]$',
			    ylabel = '$mfp_{g} \\; [m]$',
			    title = 'Ionization neutral mean free path (mid-height)')

		plot_2d(ri = ri, zi = zi, f2d = mfp_f2d, scale = "log", eq = equilibrium, title = r'Ionization mean free path')

	return (mfp_f2d, total_rate_f2d)

####################################################################################################################################################

def compute_mfp(cfg = None, verbose = False):

	location = os.path.join(cfg['baserun'],
							'input',
							cfg['run'],
							cfg['plasma']['ASTRA']['ASTRA_directory'])

	eqFile = cfg['plasma']['ASTRA']['geqdsk_name']
	eqTime = cfg['plasma']['ASTRA']['geqdsk_time']

	equilibrium = import_eqdsk(os.path.join(location, eqFile + '_' + eqTime + '.geqdsk'))

	(ng_, Tg_, np_, Tp_) = make_homemade_distributions(cfg = cfg, verbose = False)

	(ng_raw_2d, ng_f2d) = ng_
	(Tg_raw_2d, Tg_f2d) = Tg_
	(np_raw_2d, np_f2d) = np_
	(Tp_raw_2d, Tp_f2d) = Tp_

	ri = equilibrium.r_data
	zi = equilibrium.z_data

	# databases = ['AMMONX', 'AMMONX', 'AMMONX']
	# groups    = ['H.2',    'H.2',    'H.2']
	# reactions = ['R-H-H',  'R-H-H2', 'R-H2-H2']
	# vth0s     = [np.sqrt(2 * qe * Tg_raw_2d / mD), np.sqrt(2 * qe * Tg_raw_2d / mD),  np.sqrt(2 * qe * Tg_raw_2d / mD2)]
	# vth1s     = [np.sqrt(2 * qe * Tg_raw_2d / mD), np.sqrt(2 * qe * Tg_raw_2d / mD2), np.sqrt(2 * qe * Tg_raw_2d / mD2)]

	# Tg_raw_2d *= 1E-01

	vthe = np.sqrt(2 * qe * Tp_raw_2d / me)
	vthi = np.sqrt(2 * qe * Tp_raw_2d / mD)
	vthD = np.sqrt(2 * qe * Tg_raw_2d / mD)
	vthD2 = np.sqrt(2 * qe * Tg_raw_2d / mD2)

	databases = ['AMJUEL',  'AMJUEL',  'HYDHEL',  'AMJUEL',  'AMJUEL',  'AMJUEL',  'AMJUEL',  'AMJUEL',  'AMMONX',  'AMMONX',  'AMMONX']
	groups    = ['H.4',     'H.4',     'H.3',     'H.2',     'H.3',     'H.4',     'H.4',     'H.4',     'H.2',     'H.2',     'H.2']
	reactions = ['2.1.5',   '2.2.9',   '3.1.8',   '3.2.3',   '0.3T',    '2.2.5g',  '2.2.10',  '2.1.8',   'R-H-H',   'R-H-H2',  'R-H2-H2']
	n0s       = [np_raw_2d, np_raw_2d, np_raw_2d, np_raw_2d, np_raw_2d, np_raw_2d, np_raw_2d, np_raw_2d, ng_raw_2d, ng_raw_2d, ng_raw_2d]
	n1s       = [ng_raw_2d, ng_raw_2d, ng_raw_2d, ng_raw_2d, ng_raw_2d, ng_raw_2d, ng_raw_2d, np_raw_2d, ng_raw_2d, ng_raw_2d, ng_raw_2d]
	t0s       = [Tp_raw_2d, Tp_raw_2d, Tp_raw_2d, Tp_raw_2d, Tp_raw_2d, Tp_raw_2d, Tp_raw_2d, Tp_raw_2d, Tg_raw_2d, Tg_raw_2d, Tg_raw_2d]
	t1s       = [Tg_raw_2d, Tg_raw_2d, Tg_raw_2d, Tg_raw_2d, Tg_raw_2d, Tg_raw_2d, Tg_raw_2d, Tp_raw_2d, Tg_raw_2d, Tg_raw_2d, Tg_raw_2d]
	vth0s     = [vthe,      vthe,      Tg_raw_2d, vthi,      Tg_raw_2d, vthe,      vthe,      vthi,      vthD,      vthD,      vthD2]
	vth1s     = [vthD,      vthD2,     Tg_raw_2d, vthD2,     Tg_raw_2d, vthD2,     vthD2,     vthe,      vthD,      vthD2,     vthD2]

	rates = np.zeros((np_raw_2d.shape[0], np_raw_2d.shape[1], len(databases)))

	for (i, database) in enumerate(databases):

		(group, reaction, n0, n1, t0, t1, vth0, vth1) = \
		(groups[i], reactions[i], n0s[i], n1s[i], t0s[i], t1s[i], vth0s[i], vth1s[i])

		coef = read_one_reaction(database = database, group = group, reaction = reaction)

		reaction_dict = {
			'database': database,
			'group':    group,
			'reaction': reaction,
			'fit_coefficients': coef
		}
		
		(rate, mfp0, mfp1, lnsv) = \
			evaluate_rate(reaction = reaction_dict,
						  n0 = n0, t0 = t0, vth0 = vth0,
						  n1 = n1, t1 = t1, vth1 = vth1,
						  is_solps = False)
		# mfp0[np.where(mfp0 > 1E+03)] = 1E+03
		# mfp1[np.where(mfp1 > 1E+03)] = 1E+03

		rates[:,:,i] = rate

		# mfp0_f2d = Interpolate2DLinear(ri, zi, mfp0.T, extrapolate = True)
		# mfp1_f2d = Interpolate2DLinear(ri, zi, mfp1.T, extrapolate = True)

		# if cfg['semilogyting']['plot_homemade_emission'] is True and verbose is True:

		# 	plot_1d(x = ri, y = mfp0_f2d, y0 = equilibrium.magnetic_axis.y, scale = "log",
		# 		    xlabel = '$R \\; [m]$',
		# 		    ylabel = '$mfp \\; [m]$',
		# 		    title = ' - '.join(['0', database, group, reaction]))
		# 	plot_1d(x = ri, y = mfp1_f2d, y0 = equilibrium.magnetic_axis.y, scale = "log",
		# 		    xlabel = '$R \\; [m]$',
		# 		    ylabel = '$mfp \\; [m]$',
		# 		    title = ' - '.join(['1', database, group, reaction]))

		# 	plot_2d(ri = ri, zi = zi, f2d = mfp0_f2d, scale = "log", eq = equilibrium, title = ' - '.join(['0', database, group, reaction]))
		# 	plot_2d(ri = ri, zi = zi, f2d = mfp1_f2d, scale = "log", eq = equilibrium, title = ' - '.join(['1', database, group, reaction]))

	ionz = rates[:,:,0] + rates[:,:,1]
	cx   = rates[:,:,2] + rates[:,:,3]
	el   = rates[:,:,4]
	ds   = rates[:,:,5] + rates[:,:,6] 
	rc   = rates[:,:,7]
	nn   = rates[:,:,8] + rates[:,:,9] + rates[:,:,10]

	ionz_ds_probability = (ionz + ds) / (ionz + cx + el + ds + rc + nn)
	# ionz_ds_probability = (ionz) / (ionz + cx + el)
	ionz_ds_probability[np.where(ionz_ds_probability < 1E-02)] = 1E-02

	ionz_ds_probability_f2d = Interpolate2DLinear(ri, zi, ionz_ds_probability.T, extrapolate = True)

	plot_1d(x = ri, y = ionz_ds_probability_f2d, y0 = equilibrium.magnetic_axis.y, scale = "lin",
		    xlabel = '$R \\; [m]$',
		    ylabel = '$[-]$',
		    title = '(IONZ + DS) / (IONZ + CX + EL + DS + RC + NN)')
	plot_2d(ri = ri, zi = zi, f2d = ionz_ds_probability_f2d, scale = "lin", eq = equilibrium, title = '(IONZ + DS) / (IONZ + CX + EL + DS + RC + NN)')

	plt.show()

	return

####################################################################################################################################################

# def compute_mfp_average(ri = None, zi = None, mfp_f2d = None, ng_f2d = None):
def compute_mfp_average(cfg = None):

	# ri = np.linspace(6.404E-01, 7.282E-01, 100)
	ri = np.linspace(5E-01, 8E-01, 100)
	# ri = np.linspace(1.7E-01, 3E-01, 100)
	zi = np.array([-1E-02])

	(mfp_f2d, sv_f2d) = compute_mfp_ionz(cfg = cfg, verbose = False)
	((_, ng_f2d), _, (_, np_f2d), (_, Tp_f2d)) = make_homemade_distributions(cfg = cfg, verbose = False)

	# assumptions:
	# - same poloidal surface Spol of the cells
	# - 1D vertical/horizontal chord
	# - points uniformly spaced along the chord

	if ri.size > 1 and zi.size > 1:
		raise ValueError('Chord must be either vertical or horizontal.')

	if   ri.size > 1:
		is_uniform = ri == np.linspace(ri[0], ri[-1], ri.size)
		if False in is_uniform:
			raise ValueError('ri must be uniformly spaced.')
		is_horizontal = True
		size = ri.size
		zi = zi[0] * np.ones((size, 1))
	elif zi.size > 1:
		is_uniform = zi == np.linspace(zi[0], zi[-1], zi.size)
		if False in is_uniform:
			raise ValueError('zi must be uniformly spaced.')
		is_horizontal = False
		size = zi.size
		ri = ri[0] * np.ones((size, 1))

	mfp = np.zeros((size, 1))
	ng_ = np.zeros((size, 1))
	np_ = np.zeros((size, 1))
	Tp_ = np.zeros((size, 1))
	pp_ = np.zeros((size, 1))
	sv_ = np.zeros((size, 1))

	for i in range(size):
		mfp[i] = mfp_f2d(ri[i], zi[i])
		ng_[i] = ng_f2d( ri[i], zi[i])
		np_[i] = np_f2d( ri[i], zi[i])
		Tp_[i] = Tp_f2d(ri[i], zi[i])
		pp_[i] = 2 * np_[i] * Tp_[i]
		sv_[i] = sv_f2d(ri[i], zi[i])

	# integration along the 1D chord:
	#
	# - 2 * np.pi = constant
	# - Spol      = constant
	# - dl        = constant

	Ng = ri * ng_
	Np = ri * np_
	Ep = ri * pp_
	SV = ri * sv_
	mfp_average_Ng = (Ng / Ng.sum() * mfp).sum()
	mfp_average_Np = (Np / Np.sum() * mfp).sum()
	mfp_average_Ep = (Ep / Ep.sum() * mfp).sum()
	mfp_average_SV = (SV / SV.sum() * mfp).sum()

	R_average_SV  = (SV / SV.sum() * ri).sum()
	mfp_average_R = mfp[np.abs(ri - R_average_SV).argmin()]

	###############

	from raysect.core.math.random import uniform
	from statistics import stdev

	num_samples = int(1E+04)
	fp_MC = np.zeros((size, num_samples))
	mfp_MC = np.zeros((size, 1))
	mfp_SV = np.zeros((size, 1))
	err_MC = np.zeros((size, 1))
	Sigmas = np.zeros((size, 1))

	for ir in range(size):
		Sigma = sv_[ir] / np.sqrt(2 * 1.6E-19 * Tp_[ir] / 9E-31) / ng_[ir]
		Sigmas[ir] = Sigma
		for isample in range(num_samples):
			fp_MC[ir, isample] = - np.log(uniform()) / Sigma
		mfp_MC[ir] = fp_MC[ir,:].mean()

	mfp_MC[np.where(mfp_MC > 1E+03)] = 1E+03
	err_MC = mfp_MC / np.sqrt(num_samples)

	###############

	plt.figure()
	if is_horizontal: li = ri
	else:             li = zi
	plt.semilogy(li, mfp,    'k-', label = '$mfp(R)$', linewidth = 6)
	plt.semilogy(li, mfp_MC, 'k-.', label = '$mfp(R) MC$', linewidth = 6)
	plt.semilogy(li, ng_ / ng_.max() * mfp.max(),  'm-', label = 'normalised neutral density')
	plt.semilogy(li, np_ / np_.max() * mfp.max(),  'b-', label = 'normalised plasma density')
	plt.semilogy(li, pp_ / pp_.max() * mfp.max(),  'c-', label = 'normalised plasma pressure')
	plt.semilogy(li, sv_ / sv_.max() * mfp.max(),  'g-', label = 'normalised reaction rate')
	plt.semilogy([li[0], li[-1]], [mfp_average_Ng, mfp_average_Ng], 'm--', label = 'mfp neutral density average')
	plt.semilogy([li[0], li[-1]], [mfp_average_Np, mfp_average_Np], 'b--', label = 'mfp plasma density average')
	plt.semilogy([li[0], li[-1]], [mfp_average_Ep, mfp_average_Ep], 'c--', label = 'mfp plasma pressure average')
	plt.semilogy([li[0], li[-1]], [mfp_average_SV, mfp_average_SV], 'g--', label = 'mfp reaction rate average')
	plt.semilogy([li[0], li[-1]], [mfp_average_R, mfp_average_R], 'k--', label = 'mfp radius average')
	plt.xlabel('$R \\; [m]$')
	plt.ylabel('$[-]$')
	plt.legend(loc = 'center left')
	plt.ylim(bottom = 1E-02, top = 2E+03)
	plt.show()

	return

