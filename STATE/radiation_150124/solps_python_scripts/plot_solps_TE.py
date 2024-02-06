import os
import time
import pickle
import numpy as np
import matplotlib as mtplb
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

from decimal import Decimal

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from solps_python_scripts.utilities.units import units
from solps_python_scripts.utilities.names import name
from solps_python_scripts.utilities.last10s import read_last10s
from solps_python_scripts.utilities.plotting import make_ticks, triplot, to_triangles
from solps_python_scripts.utilities.accessories import load_pickle, find_sp44, populate_b2
from solps_python_scripts.utilities.accessories import rearrange_quadrangles, rearrange_triangles

from solps_python_scripts.read_b2fgmtry import read_b2fgmtry
from solps_python_scripts.read_b2fstate import read_b2fstate
from solps_python_scripts.read_b2fplasmf import read_b2fplasmf
from solps_python_scripts.read_ft44 import read_ft44
from solps_python_scripts.read_ft46 import read_ft46
from solps_python_scripts.read_triangle_mesh import read_triangle_mesh

from solps_python_scripts.reactions.compute_rates import compute_rates

from st40_phys_viewer import Get_data

#########################################################################################################

def divide_by_vol(values = None, vol = None):

	# Author: Matteo Moscheni
	# E-mail: matteo.moscheni@tokamakenergy.co.uk
	# February 2022

	try:		
		for i in range(values.shape[2]):
			values[:,:,i] /= vol
	except:
		values /= vol

	return values

#########################################################################################################

def divide_by_e(values = None):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	return values / 1.6E-19

#########################################################################################################

def plot_vs_ASTRA(where = ".", what = None,
		  		  shot_number = None, ASTRA_run = None, time = None,
			   	  vmin = None, vmax = None, scale = "lin",
	           	  legend = True, labels = None, loc = "best", grid = True, save = False):
	"""

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # June 2022

    Input example:
    >>> shot_number = 13110014
    >>> ASTRA_run   = 'ASTRA#RUN700'

    """

	if not isinstance(where, list):
		where = [where]
	if not isinstance(what, list):
		what = [what]
	if labels is not None and not isinstance(labels, list):
		labels = [labels]

	#####
	# get ASTRA data from MDSplus tree

	ASTRA = Get_data(shot_number, ASTRA_run)

	if isinstance(ASTRA.get('TIME'), str) and ASTRA.get('TIME') == 'FAILED':
		raise ValueError('Shot number {} - {} - does NOT exist!'.format(shot_number, ASTRA_run))

	time_ASTRA = ASTRA.get('TIME')
	time_index = int(np.where(np.abs(time_ASTRA - time) == np.abs(time_ASTRA - time).min())[0][0])

	# select OMP half of midplane coordinate at specified time
	# and shift it to get (R - Rsep) ~ SOLPS convention
	R_ASTRA = ASTRA.get('R_MID_PROFILES.R')
	R_ASTRA = R_ASTRA[time_index, int(R_ASTRA.shape[1] / 2):] - R_ASTRA[time_index,-1]

	#####

	last10s = {}
	marker = ['o', 's', 'x', '^', 'v', '*']

	for this in range(len(what)):

		for i in range(len(where)):

			here = where[i]

			last10s[here] = {}
			last10s[here] = read_last10s(where = here, save = True)
			# try:    b2fstate = load_pickle(where = here, what = "b2fstate")
			# except: b2fstate = read_b2fstate(where = here, save = True)
			b2fstate = read_b2fstate(where = here, save = True)

		# variables allowed (both OMP and IMP)
		#
		# - ne
		# - ni
		# - te
		# - ti
		# - chii
		# - chie
		# - diff
		# - Qi vs. Psol
		# - Qe vs. Psol
		# - Q_RAD vs. total radiated SOLPS (see CHERAB)

		variable_SOLPS = what[this]

		if variable_SOLPS[-1] != 'a' and variable_SOLPS[-1] != 'i':
			raise ValueError('Comparing ASTRA with non-midplane profiles is not allowed...')
		elif '3d' in variable_SOLPS:
			# usual convention
			if variable_SOLPS[:-1] == 'na3d':
				raise ValueError('Multiple species not implemented yet...')
			variable_ASTRA = variable_SOLPS[:2].upper()
		else:
			# transport coefficient convention
			variable_ASTRA = variable_SOLPS.upper()

		try:
			# get ASTRA data for specified timestep
			data_ASTRA = ASTRA.get('PROFILES.' + variable_ASTRA)[time_index,:]
		except:
			raise ValueError('{} not found in ASTRA results...'.format(variable_ASTRA))

		plt.figure()

		for (i, here) in enumerate(where):
			
			data = last10s[here][what[this]]

			if   units(of = what[this]) == 'W \\rightarrow W \\cdot m^{-3}':
				vol = last10s[here]['vol3' + what[this][-2:]]
				data[:,1] = divide_by_vol(values = data[:,1], vol = vol[:,1])
			elif units(of = what[this]) == 'A \\rightarrow m^{-3}':
				data = divide_by_e(values = data)

			if i == 0:
				colors = pylab.cm.gnuplot2(np.linspace(0, 1, (data.shape[1]-1) * len(where) + 1))
			
			# if data.shape[1] > 1:
			species = b2fstate['species']
			for j in range(1,data.shape[1]):
				jcol = (j-1) + i * (data.shape[1]-1)
				if scale == "lin":
					if i == 0 and data.shape[1] > 2:
						plt.plot(data[:,0], data[:,j], marker[i] + '-', color = colors[jcol], label = species[j-1])
						plt.plot(R_ASTRA, data_ASTRA, marker[i] + 'b-')
					else:
						if labels is not None:
							plt.plot(data[:,0], data[:,j], marker[i] + '-', color = colors[jcol], label = labels[i]) 
						else:
							plt.plot(data[:,0], data[:,j], marker[i] + '-', color = colors[jcol])     
							plt.plot(R_ASTRA, data_ASTRA, marker[i] + 'b-')
				elif scale == "log":
					if i == 0 and data.shape[1] > 2:
						plt.semilogy(data[:,0], data[:,j], marker[i] + '-', color = colors[jcol], label = species[j-1])
						plt.semilogy(R_ASTRA, data_ASTRA, marker[i] + 'b-')
					else:
						if labels is not None:
							plt.semilogy(data[:,0], data[:,j], marker[i] + '-', color = colors[jcol], label = labels[i])
						else:
							plt.semilogy(data[:,0], data[:,j], marker[i] + '-', color = colors[jcol])
							plt.semilogy(R_ASTRA, data_ASTRA, marker[i] + 'b-')
		
		plt.xlabel(r'$s \; [m]$')
		plt.ylabel(r"$[" + units(of = what[this]) + "]$")
		if vmin is not None: plt.ylim(bottom = vmin)
		if vmax is not None: plt.ylim(top = vmax)
		plt.title(name(of = what[this]) or what[this])
		if legend is True: plt.legend(loc = loc)
		plt.grid(grid)
	
	plt.show()

	return

#########################################################################################################