import os
import time
import pickle
import numpy as np

import matplotlib as mtplb
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from solps_python_scripts.read_b2mn import read_b2mn
from solps_python_scripts.read_b2fgmtry import read_b2fgmtry
from solps_python_scripts.utilities.accessories import load_pickle

from solps_python_scripts.utilities.units import units
from solps_python_scripts.utilities.names import name

#########################################################################################################

def make_ticks(vmin_log10 = None, vmax_log10 = None):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	num_ticks = 6

	ticks = np.logspace(vmin_log10, vmax_log10, num_ticks)
	tick_list = []
	tick_labels = []

	for tick in ticks:
		tick_list   += [tick]
		tick_labels += ["{:.1e}".format(tick)]

	return (np.log10(np.array(tick_list)), tick_labels)

#########################################################################################################

def plot_separatrix(where = '.', colour = 'w', annotate_ix = False):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	# try:    b2fgmtry = load_pickle(where = where, verbose = False, what = "b2fgmtry")
	# except: b2fgmtry = read_b2fgmtry(where = where, verbose = False, save = True)
	b2fgmtry = read_b2fgmtry(where = where, verbose = False, save = True)

	iy = int(b2fgmtry['ny'] / 2)

	if len(b2fgmtry['rightcut']) == 2:

		ix_mid = int((b2fgmtry['rightcut'][1] + b2fgmtry['leftcut'][1]) / 2 - 1)

		for ix in range(ix_mid - 2):

			x01 = [b2fgmtry['crx'][ix,iy,0], b2fgmtry['crx'][ix,iy,1]]
			y01 = [b2fgmtry['cry'][ix,iy,0], b2fgmtry['cry'][ix,iy,1]]
			plt.plot(x01, y01, colour + '-')
			if annotate_ix is True:
				plt.annotate(str(ix), (np.mean(x01), np.mean(y01)))

		for ix in range(ix_mid, b2fgmtry['nx']):

			x01 = [b2fgmtry['crx'][ix,iy,0], b2fgmtry['crx'][ix,iy,1]]
			y01 = [b2fgmtry['cry'][ix,iy,0], b2fgmtry['cry'][ix,iy,1]]
			plt.plot(x01, y01, colour + '-')
			if annotate_ix is True:
				plt.annotate(str(ix), (np.mean(x01), np.mean(y01)))

	else:

		for ix in range(b2fgmtry['nx']):

			x01 = [b2fgmtry['crx'][ix,iy,0], b2fgmtry['crx'][ix,iy,1]]
			y01 = [b2fgmtry['cry'][ix,iy,0], b2fgmtry['cry'][ix,iy,1]]
			plt.plot(x01, y01, colour + '-')
			if annotate_ix is True:
				plt.annotate(str(ix), (np.mean(x01), np.mean(y01)))

	return

#########################################################################################################

def triplot(where = None, what = None, value = None, sp = None, nodes = None, cells = None, vmin = None, vmax = None,
	        cmap = None, scale = "log", top = None, bottom = None, right = None, left = None, mask = True, save = False):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	fluid_neutral_not_plotted = False
	zeroes_not_plotted         = False

	# skip neutral fluid species if EIRENE is used
	b2mn = read_b2mn(where = where)

	for i in range(len(sp)):

		try:    values_original = value[:,i]
		except: values_original = value

		if (values_original > 0).max() == True:

			zeroes_not_plotted = True

			if values_original[np.isnan(values_original) == 0].min() < 0:
				signs = [+1, -1]
			else:
				signs = [+1]

			for sign in signs:

				index_mask = values_original == np.nan

				values = sign * values_original
				values[values < 0] = 0

				index_mask += values == 0

				if vmin is not None:
					# values[values < vmin] = vmin
					if scale == "log": vmin_log10 = np.log10(vmin)
					else: vmin_log10 = vmin
				if vmax is not None:
					values[values > vmax] = vmax
					if scale == "log": vmax_log10 = np.log10(vmax)
					else: vmax_log10 = vmax
				
				if scale == "log":
					values = np.log10(values + 1E-15)
					if vmin is None:
						vmin_log10 = values[values > -14].min()
				elif vmin is None: vmin_log10 = values.min()
				if vmax is None: vmax_log10 = values.max()

				# plot neutrals (0+) only if EIRENE = OFF

				if b2mn['b2mndr_eirene'] == 1 and '0+' in sp[i]: pass

				else:

					fluid_neutral_not_plotted = True
					fig, ax = plt.subplots()

					triangulation = tri.Triangulation(nodes[:,0], nodes[:,1], cells)
					if mask is True: triangulation.set_mask(index_mask)
					
					tpc = ax.tripcolor(triangulation, values,
					 				   shading = 'flat', cmap = cmap or 'jet',
					 				   vmin = vmin_log10, vmax = vmax_log10)

					plot_separatrix(where = where)

					cbar = fig.colorbar(tpc, label = '$[' + units(of = what) + ']$')
					(ticks, tick_labels) = make_ticks(vmin_log10, vmax_log10)
					if scale == "log":
						cbar.set_ticks(ticks)
						cbar.set_ticklabels(tick_labels)
					if top is not None: ax.set_ylim(top = top)
					if bottom is not None: ax.set_ylim(bottom = bottom)
					if left is not None: ax.set_xlim(left = left)
					if right is not None: ax.set_xlim(right = right)
					ax.set_xlabel(r'$r\;[m]$')		
					ax.set_ylabel(r'$z\;[m]$')
					if sign == +1:
						try: ax.set_title(sp[i] + " " + name(of = what) or sp[i] + " " + what)
						except: ax.set_title(what)
					elif sign == -1:
						try: ax.set_title("$(" + str(sign) + ") \\cdot$ " + sp[i] + " " + name(of = what) or sp[i] + " " + what)
						except: ax.set_title("$(" + str(sign) + ") \\cdot$ " + what)
					ax.set_aspect('equal')

					if save is True: plt.savefig(os.path.join(where, what + '_' + str(sp[i]) + '_' + str(sign) + '.png'))
		else:
			pass

	if zeroes_not_plotted is False:
		print()
		print('ATTENTION: nothing to plot here, all zeroes!')
	elif fluid_neutral_not_plotted is False:
		print()
		print('ATTENTION: nothing to plot here, all fluid neutrals but EIRENE used!')

	return

#########################################################################################################

def to_triangles(nsp = None, vs = None, b2fgmtry = None):

	"""

	to_triangles create a triangulation starting from the quadrangular
	B2 mesh (to simplify plotting and homogenise with plot_eirene)

	From each B2 quadrangle, 2 triangles are obtained (EIRENE-mesh-like)
	=> each value in the quadrangular cell centre must be duplicated

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	"""

	cells = []
	nodes = []
	value = []
	inode = 0

	for ix in range(b2fgmtry['nx']):
		for iy in range(b2fgmtry['ny']):

			# first triangle

			nodes += [[b2fgmtry['crx'][ix,iy,1], b2fgmtry['cry'][ix,iy,1]]]
			nodes += [[b2fgmtry['crx'][ix,iy,3], b2fgmtry['cry'][ix,iy,3]]]
			nodes += [[b2fgmtry['crx'][ix,iy,2], b2fgmtry['cry'][ix,iy,2]]]

			cells += [[inode, inode + 1, inode + 2]]

			inode += 3

			# second triangle

			nodes += [[b2fgmtry['crx'][ix,iy,2], b2fgmtry['cry'][ix,iy,2]]]
			nodes += [[b2fgmtry['crx'][ix,iy,0], b2fgmtry['cry'][ix,iy,0]]]
			nodes += [[b2fgmtry['crx'][ix,iy,1], b2fgmtry['cry'][ix,iy,1]]]

			cells += [[inode, inode + 1, inode + 2]]

			inode += 3

			# duplicate values

			if nsp > 1:
				add = np.zeros(nsp)
				for isp in range(vs.shape[2]):
					add[isp] = vs[ix,iy,isp]
				value += [add]
				value += [add] # 2 triangles in 1 quadrangle
			else:
				value += [vs[ix,iy]]
				value += [vs[ix,iy]] # 2 triangles in 1 quadrangle

	return (np.array(value), np.array(nodes), np.array(cells))
