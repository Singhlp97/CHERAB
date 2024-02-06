import os
import csv
import numpy as np
import matplotlib.pyplot as plt

from cherab.tools.equilibrium import import_eqdsk
from cherab.core.math import Interpolate1DLinear, Interpolate2DLinear

from import_SOLPS_data import load_solps_data
from astra_simulation import ASTRASimulation

##############################################################################################
##############################################################################################
##############################################################################################

def load_astra_data(cfg = None):

	location = os.path.join(cfg['baserun'],
							'input',
							cfg['run'],
							cfg['plasma']['ASTRA']['ASTRA_directory'])

	ASTRAfile = cfg['plasma']['ASTRA']['ASTRA_file']
	eqFile    = cfg['plasma']['ASTRA']['geqdsk_name']
	eqTime    = cfg['plasma']['ASTRA']['geqdsk_time']

	# MMM: different weights
	use_weighted_emission = cfg['plasma']['SOLPS']['weights']['use_weighted_emission']
	if use_weighted_emission is True:
		# only weight of main emission by assumption (ASTRA from main emission)
		weight = cfg['plasma']['SOLPS']['weights']['weight_main_emission']

	# MMM R-dependent spatial rescaling to match SOLPS
	use_ASTRA_rescaling = cfg['plasma']['ASTRA']['use_ASTRA_rescaling']

	plot_ASTRA_emission = cfg['plotting']['plot_ASTRA_emission'] 

	##############################################################################################
	##############################################################################################
	##############################################################################################

	ASTRA_simulation = ASTRASimulation()

	equilibrium = import_eqdsk(os.path.join(location, eqFile + "_" + eqTime + ".geqdsk"))
	ASTRA_simulation.equilibrium = equilibrium

	r = equilibrium.r_data
	z = equilibrium.z_data

	data = np.loadtxt(open(location + ASTRAfile, "r"), delimiter = ",", skiprows = 1)

	# properly rescaling
	# if data[inside_lcfs,3].max() < 1E+10:
	# 	scale = 1E+14
	# else:
	# 	scale = 1.0
	scale = 1.0

	psiN = data[:,2]
	inside_lcfs = (psiN <= 1.0)
	psiN = psiN[inside_lcfs]

	a = data[inside_lcfs,1]                      # [m]
	ne = data[inside_lcfs,3] * 1E+19             # [1 / m^3]
	Te = data[inside_lcfs,4] * 1E+03             # [eV]
	n0 = data[inside_lcfs,5] * 1E+19            # [1 / m^3]
	Ha = data[inside_lcfs,6] * scale * 1E+06 # [photons / s / m^3]

	# extrapolate = True => need to manually set to zero outside range
	neVsPsiN = Interpolate1DLinear(psiN, ne, extrapolate = True)
	TeVsPsiN = Interpolate1DLinear(psiN, Te, extrapolate = True)
	n0VsPsiN = Interpolate1DLinear(psiN, n0, extrapolate = True)
	HaVsPsiN = Interpolate1DLinear(psiN, Ha, extrapolate = True)

	# (r,z) from eqdsk => psi(r,z) from eqdsk => psiN(psi(r,z)) from formula => Ha(psiN(psi(r,z))) from ASTRA
	#
	# formula: pfi_norm = (psi - psi_axis) / (psi_lcfs - psi_axis)

	def normalise_psi(eq = None):
		return (eq.psi_data - eq.psi_axis) / (eq.psi_lcfs - eq.psi_axis)

	psi_normalised_data = normalise_psi(eq = equilibrium)

	##############################################################################################
	##############################################################################################
	##############################################################################################

	if use_ASTRA_rescaling is True:
		L  = cfg['plasma']['ASTRA']['ASTRA_rescaling_L']
		r0 = cfg['plasma']['ASTRA']['ASTRA_rescaling_r0']
	else:
		L = np.inf
		r0 = 0.0

	def ASTRA_rescale(r = None, r0 = r0, L = L):
		return np.exp(-(r - r0) / L)

	# mind index convention
	ne = np.zeros((np.size(z), np.size(r)))
	Te = np.zeros((np.size(z), np.size(r)))
	n0 = np.zeros((np.size(z), np.size(r)))
	Ha = np.zeros((np.size(z), np.size(r)))
	for i in range(Ha.shape[0]):
		for j in range(Ha.shape[1]):
			# emission = 0.0 outside lcfs
			inside = equilibrium.inside_lcfs(r[j], z[i])
			ne[i,j] = neVsPsiN(psi_normalised_data[j,i]) * inside
			Te[i,j] = TeVsPsiN(psi_normalised_data[j,i]) * inside
			n0[i,j] = n0VsPsiN(psi_normalised_data[j,i]) * inside
			Ha[i,j] = HaVsPsiN(psi_normalised_data[j,i]) * inside * ASTRA_rescale(r = r[j], r0 = r0, L = L)

	ASTRA_simulation.electron_density_raw = ne
	ASTRA_simulation.electron_temperature_raw = Te
	ASTRA_simulation.neutral_density_raw = n0
	ASTRA_simulation.Ha_emission_raw = Ha

	if use_weighted_emission is True: ASTRA_simulation.Ha_emission_raw *= weight

	################

	SOLPSsim = load_solps_data(cfg = cfg)

	################

	# CAUTION pt. 1
	#
	# Interpolate2DCubic would MISTAKENLY give:
	# - overestimated maximum
	# - underestimated minimum < 0.0
	#
	# => Interpolate2Dlinear must be used!!
	#
	# CAUTION pt. 2
	#
	# Assumption: entire SOLPS-ITER result, ASTRA's cut where SOLPS-ITER's end in the core

	where_non_zero = (1.0 - SOLPSsim._inside_mesh) * equilibrium.inside_lcfs

	ne_f2d = Interpolate2DLinear(r, z, ne.T, extrapolate = True) * where_non_zero
	Te_f2d = Interpolate2DLinear(r, z, Te.T, extrapolate = True) * where_non_zero
	n0_f2d = Interpolate2DLinear(r, z, n0.T, extrapolate = True) * where_non_zero
	Ha_f2d = Interpolate2DLinear(r, z, Ha.T, extrapolate = True) * where_non_zero

	ASTRA_simulation.electron_density_f2d = ne_f2d
	ASTRA_simulation.electron_temperature_f2d = Te_f2d
	ASTRA_simulation.neutral_density_f2d = n0_f2d
	ASTRA_simulation.Ha_emission_f2d = Ha_f2d

	##############################################################################################
	##############################################################################################
	##############################################################################################

	if plot_ASTRA_emission is True:

		##############################################################################################

		plt.figure()
		plt.plot(a, data[inside_lcfs,6] * scale * 1E+06, 'ko-')
		plt.xlabel('$r\:[m]$')
		plt.ylabel(r'$H_{\alpha}\;[ph/s/m^3]$')

		if cfg['plotting']['save_figures'] is True:
			plt.savefig(os.path.join(cfg['output_directory_extended'], 'ASTRA_Ha_1D.png'))

		##############################################################################################

		fig, ax = plt.subplots()
		c = ax.pcolormesh(r, z, Ha)
		fig.colorbar(c, ax = ax)
		plt.plot(equilibrium.lcfs_polygon[:,0], equilibrium.lcfs_polygon[:,1], 'm-')
		ax.set_title(r'ASTRA: $H_{\alpha}\;[ph\cdot s^{-1}\cdot m^{-3}]$ - original')
		ax.axis('equal')
		ax.axis([0.1, 0.7, -0.35, 0.35])

		if cfg['plotting']['save_figures'] is True:
			plt.savefig(os.path.join(cfg['output_directory_extended'], 'ASTRA_Ha_complete.png'))

		##############################################################################################

		num = 250
		ri = np.linspace(r.min(), r.max(), num)
		zi = np.linspace(z.min(), z.max(), num)

		# mind imshow convention
		ne = np.zeros((np.size(zi), np.size(ri)))
		Te = np.zeros((np.size(zi), np.size(ri)))
		n0 = np.zeros((np.size(zi), np.size(ri)))
		Ha = np.zeros((np.size(zi), np.size(ri)))
		Ha_check = np.zeros((np.size(ri))) # compute line-integral value to compare against line diodes
		for i in range(Ha.shape[0]):
			for j in range(Ha.shape[1]):
				# emission = 0.0 outside lcfs
				ne[i,j] = ne_f2d(ri[j], zi[i])
				Te[i,j] = Te_f2d(ri[j], zi[i])
				n0[i,j] = n0_f2d(ri[j], zi[i])
				Ha[i,j] = Ha_f2d(ri[j], zi[i])
				if i == int(np.round(Ha.shape[0] / 2)):
					Ha_check[j] = Ha_f2d(ri[j], zi[i])

		Ha_integral = ((ri[1:] - ri[:-1]) * (Ha_check[1:] + Ha_check[:-1]) * 0.5).sum()
		# print(Ha_integral, Ha_integral / 4 / np.pi)
		# plt.figure()
		# plt.plot(ri, Ha_check)
		# plt.show()

		##############################################################################################

		fig, ax = plt.subplots()
		c = ax.pcolormesh(ri, zi, ne)
		fig.colorbar(c, ax = ax)
		plt.plot(equilibrium.lcfs_polygon[:,0], equilibrium.lcfs_polygon[:,1], 'm-')
		ax.set_title(r'ASTRA: $n_e\;[m^{-3}]$')
		ax.axis('equal')
		ax.axis([0.1, 0.7, -0.35, 0.35])

		if cfg['plotting']['save_figures'] is True:
			plt.savefig(os.path.join(cfg['output_directory_extended'], 'ASTRA_ne.png'))

		##############################################################################################

		fig, ax = plt.subplots()
		c = ax.pcolormesh(ri, zi, Te)
		fig.colorbar(c, ax = ax)
		plt.plot(equilibrium.lcfs_polygon[:,0], equilibrium.lcfs_polygon[:,1], 'm-')
		ax.set_title(r'ASTRA: $T_e\;[eV]$')
		ax.axis('equal')
		ax.axis([0.1, 0.7, -0.35, 0.35])

		if cfg['plotting']['save_figures'] is True:
			plt.savefig(os.path.join(cfg['output_directory_extended'], 'ASTRA_Te.png'))

		##############################################################################################

		fig, ax = plt.subplots()
		c = ax.pcolormesh(ri, zi, n0)
		fig.colorbar(c, ax = ax)
		plt.plot(equilibrium.lcfs_polygon[:,0], equilibrium.lcfs_polygon[:,1], 'm-')
		ax.set_title(r'ASTRA: $n_{H_0}\;[m^{-3}]$')
		ax.axis('equal')
		ax.axis([0.1, 0.7, -0.35, 0.35])

		if cfg['plotting']['save_figures'] is True:
			plt.savefig(os.path.join(cfg['output_directory_extended'], 'ASTRA_nH.png'))

		##############################################################################################

		fig, ax = plt.subplots()
		c = ax.pcolormesh(ri, zi, Ha)
		fig.colorbar(c, ax = ax)
		# plt.plot(equilibrium.lcfs_polygon[:,0], equilibrium.lcfs_polygon[:,1], 'm-')
		ax.set_title(r'ASTRA: $H_{\alpha}\;[ph\cdot s^{-1}\cdot m^{-3}]$')
		ax.axis('equal')
		ax.axis([0.1, 0.7, -0.35, 0.35])

		if cfg['plotting']['save_figures'] is True:
			plt.savefig(os.path.join(cfg['output_directory_extended'], 'ASTRA_Ha_2D.png'))

		##############################################################################################

		#plt.show()

	##############################################################################################

	return ASTRA_simulation