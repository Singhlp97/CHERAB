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
from solps_python_scripts.utilities.accessories import load_pickle, find_sp44, populate_b2, compute_integral
from solps_python_scripts.utilities.accessories import rearrange_quadrangles, rearrange_triangles

from solps_python_scripts.read_b2mn import read_b2mn
from solps_python_scripts.read_b2fgmtry import read_b2fgmtry
from solps_python_scripts.read_b2fstate import read_b2fstate
from solps_python_scripts.read_b2fplasmf import read_b2fplasmf
from solps_python_scripts.read_ft44 import read_ft44
from solps_python_scripts.read_ft46 import read_ft46
from solps_python_scripts.read_triangle_mesh import read_triangle_mesh

from solps_python_scripts.reactions.compute_rates import compute_rates

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

def compute_variation(values = None, is_1d = None):

	# Author: Matteo Moscheni
	# E-mail: matteo.moscheni@tokamakenergy.co.uk
	# February 2022

	# format:
	# - first column: abscissa
	# - second column onwards: different data

	if is_1d is True: index = 1
	else:             index = 0

	variation = np.zeros((values.shape[0], index + 1))
	if is_1d is True: variation[:,0] = values[:,0]

	for i in range(values.shape[0]):
		variation[i,index]  = np.max(values[i,index:]) - np.min(values[i,index:])
		if np.abs(values[i,index:].mean()) > 0:
			variation[i,index] /= np.abs(values[i,index:].mean())
		else:
			variation[i,index] = 0.0

	return variation / (values.shape[1] - index)

#########################################################################################################

def plot_1d(where = ".", what = None, vmin = None, vmax = None, scale = "lin", integral = False,
	        legend = False, labels = None, loc = "best", grid = True, variation = False, save = False):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	if not isinstance(where, list):
		if variation is True:
			raise ValueError('where must be a list if variation is True.')
		else:
			where = [where]
	if not isinstance(what, list):
		what = [what]
	if labels is not None and not isinstance(labels, list):
		labels = [labels]
	if integral is True: legend = True

	last10s = {}
	marker = ['o', 's', 'x', '^', 'v', '*']

	for this in range(len(what)):

		for i in range(len(where)):

			here = where[i]

			last10s[here] = {}
			last10s[here] = read_last10s(where = here, save = True)
			b2fstate = read_b2fstate(where = here, save = True)

		plt.figure()

		for (i, here) in enumerate(where):
			
			data = last10s[here][what[this]]

			if variation is True:
				if data.shape[1] > 2:
					raise ValueError('variation for multiple species not implemented yet.')
				if i == 0:
					cases = np.zeros((data.shape[0], len(where) + 1))
					cases[:,0] = data[:,0]
				cases[:,i+1] = data[:,1]

			if   units(of = what[this]) == 'W \\rightarrow W \\cdot m^{-3}':
				vol = last10s[here]['vol3' + what[this][-2:]]
				data[:,1] = divide_by_vol(values = data[:,1], vol = vol[:,1])
			elif units(of = what[this]) == 'A \\rightarrow m^{-3}':
				data[:,1] = divide_by_e(values = data[:,1])

			if i == 0:
				colors = pylab.cm.gnuplot2(np.linspace(0, 1, (data.shape[1]-1) * len(where) + 1))

			if variation is False:

				# data[:,1:] *= 1e-6

				# if data.shape[1] > 1:
				species = b2fstate['species']
				for j in range(1,data.shape[1]):
					jcol = (j-1) + i * (data.shape[1]-1)
					if scale == "lin":
						if i == 0 and data.shape[1] > 2:
							plt.plot(data[:,0], data[:,j], marker[i] + '-', color = colors[jcol], label = species[j-1])
						else:
							if legend is True:
								try:    label = labels[i]
								except: label = ''
								if integral is True:
									b2fgmtry = read_b2fgmtry(where = here, verbose = False)
									integ = compute_integral(b2fgmtry = b2fgmtry, last10s = last10s[here], what = what[this])
									label += ' (' + '%.3E' % Decimal(str(integ)) + ')'
								plt.plot(data[:,0], data[:,j], marker[i] + '-', color = colors[jcol], label = label)
							else:
								plt.plot(data[:,0], data[:,j], marker[i] + '-', color = colors[jcol])     
					elif scale == "log":
						if i == 0 and data.shape[1] > 2:
							plt.semilogy(data[:,0], data[:,j], marker[i] + '-', color = colors[jcol], label = species[j-1])
						else:
							if legend is True:
								try:    label = labels[i]
								except: label = ''
								if integral is True:
									b2fgmtry = read_b2fgmtry(where = here, verbose = False)
									integ = compute_integral(b2fgmtry = b2fgmtry, last10s = last10s[here], what = what[this])
									label += ' (' + '%.3E' % Decimal(str(integ)) + ')'
								plt.semilogy(data[:,0], data[:,j], marker[i] + '-', color = colors[jcol], label = label)
							else:
								plt.semilogy(data[:,0], data[:,j], marker[i] + '-', color = colors[jcol])

		if variation is True:
			data = compute_variation(values = cases, is_1d = True)
			what[this] += " relative variation"
			if scale == "lin":
				plt.plot(data[:,0], data[:,1], 'o-', color = colors[0])     
			elif scale == "log":
				plt.semilogy(data[:,0], data[:,1], 'o-', color = colors[0])

		plt.xlabel(r'$(y-y_{sep}) \; [m]$')
		plt.ylabel(r"$[" + units(of = what[this]) + "]$")
		if vmin is not None: plt.ylim(bottom = vmin)
		if vmax is not None: plt.ylim(top = vmax)
		plt.title(name(of = what[this]) or what[this])
		if legend is True: plt.legend(loc = loc)
		plt.grid(grid)
	
	plt.show()

	return

#########################################################################################################

def plot_wall_loads(where = ".", what = "wldnek", vmin = None, vmax = None, scale = "lin",
	        		legend = True, labels = None, loc = "best", grid = True, save = False):
	
	(neut, wld) = read_ft44(where = where, save = True)

	# user to make sure the plotting dimensions are coheret (see SOLPS-ITER manual)

	poly = wld['poly']
	area = wld['wlarea']
	data = wld[what] + wld['wldnep']

	fig, ax = plt.subplots(1,2)

	for i in range(poly.shape[1]):
		ax[0].plot([poly[0,i], poly[2,i]], [poly[1,i], poly[3,i]], 'ko-')
		ax[0].annotate(str(i), (poly[0,i], poly[1,i]),
			        textcoords = "offset points",
                    xytext = (0,8),
                    ha = 'left')

	ax[0].set_xlabel('$R \\; [m]$')
	ax[0].set_ylabel('$Z \\; [m]$')
	ax[0].set_title('Standard surfaces')
	ax[0].set_aspect('equal')

	ax[1].plot(np.linspace(0, poly.shape[1]-1, poly.shape[1]), data[:poly.shape[1],0] / area[:poly.shape[1]] * 1E-03, 'r-')
	ax[1].set_xlabel('element number [-]')
	ax[1].set_ylabel('$[kW \\cdot m^{-2}]$')
	ax[1].set_title(what)

	plt.show()

	return

#########################################################################################################

def plot_2d(where = ".", what = None, vmin = None, vmax = None, cmap = None, scale = "log",
	        top = None, bottom = None, right = None, left = None, mask = True, variation = False, save = False,
	        xgc1 = False, is_usn = False):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	if not isinstance(where, list):
		if variation is True:
			raise ValueError('where must be a list if variation is True.')
		else:
			where = [where]
	if not isinstance(what, list):  what  = [what]
	if not isinstance(vmin, list):  vmin  = [vmin]
	if not isinstance(vmax, list):  vmax  = [vmax]

	for j in range(len(what)):

		for i in range(len(where)):

			here = where[i]
			try: this = what[j]
			except: this = what[0]
			try: low = vmin[j]
			except: low = vmin[0]
			try: up = vmax[j]
			except: up = vmax[0]

			b2mn = read_b2mn(where = here)
			b2fgmtry = read_b2fgmtry(where = here, save = True)
			b2fstate = read_b2fstate(where = here, save = True)

			try:
				b2fplasmf = read_b2fplasmf(where = here, save = True)
			except:
				print()
				print('NO b2fplasmf > b2fstate used instead')
				print()
				b2fplasmf = b2fstate

			fort44 = {}
			fort46 = {}

			if b2mn['b2mndr_eirene'] == 1: 

				[neut, wld] = read_ft44(where = here, save = True)
				fort44['neut'] = neut
				fort44['wld']  = wld

				fort46 = read_ft46(where = here, save = True)
			
			try: 
				(sp, value) = populate_b2(what = this, b2mn = b2mn, b2fstate = b2fstate, b2fplasmf = b2fplasmf, fort44 = fort44)
				is_b2 = True
			except:
				try: value = fort46[this]
				except: raise ValueError(this + ' does NOT exist :(')
				sp = find_sp44(what = this, fort44 = fort44)
				is_b2 = False

			if is_usn is True and is_b2 is True: value = value[::-1,:]

			try:    nsp = value.shape[2]
			except: nsp = 1

			if   units(of = this) == 'W \\rightarrow W \\cdot m^{-3}':
				vol = b2fgmtry['vol']
				value = divide_by_vol(values = value, vol = vol)
			elif units(of = this) == 'A \\rightarrow m^{-3}':
				value = divide_by_e(values = value)

			if is_b2 is True:
				(value, nodes, cells) = to_triangles(nsp = nsp, vs = value, b2fgmtry = b2fgmtry)
			else:
				triangles = read_triangle_mesh(where = here, verbose = False, save = False)
				cells = triangles['cells']
				nodes = triangles['nodes']

			if variation is True:
				if i == 0: cases = np.zeros((value.shape[0], len(where)))
				try:
					if value.shape[1] > 1:
						raise ValueError('variation for multiple species not implemented yet.')
					cases[:,i] = value[:,0]
				except IndexError:
					cases[:,i] = value

			else:

				if this == "ti" or "tot" in this: sp = [""]

				triplot(where = here, what = this, value = value, sp = sp, nodes = nodes, cells = cells,
					    vmin = low, vmax = up, cmap = cmap, scale = scale, top = top, bottom = bottom,
					    right = right, left = left, mask = mask, save = save)

		if variation is True:

			value = compute_variation(values = cases, is_1d = False)
			this += ' relative variation'

			triplot(where = here, what = this, value = value, sp = sp, nodes = nodes, cells = cells,
				    vmin = low, vmax = up, cmap = cmap, scale = scale, top = top, bottom = bottom,
				    right = right, left = left, mask = mask, save = save)

	plt.show()

	return

#########################################################################################################

def plot_mesh(where = ".", what = ["b2", "eirene"],
	   		  top = None, bottom = None, right = None, left = None, save = False):

	if isinstance(what, list) is False: what = [what]

	for this in what:

		if this == "b2":

			b2fgmtry = read_b2fgmtry(where = where, verbose = False)
			
			nx  = b2fgmtry['nx']
			ny  = b2fgmtry['ny']
			crx = b2fgmtry['crx']
			cry = b2fgmtry['cry']

			(crx, cry) = rearrange_quadrangles(crx = crx, cry = cry)

			fig, ax = plt.subplots()
			for ix in range(nx):
				for iy in range(ny):
					ax.plot(crx[ix,iy,:], cry[ix,iy,:], 'b-', linewidth = 0.5)

		elif this == "eirene":

			triangles = read_triangle_mesh(where = where, verbose = False, save = False)
			cells = triangles['cells']
			nodes = triangles['nodes']

			cells = rearrange_triangles(cells = cells)

			fig, ax = plt.subplots()
			for ic in range(cells.shape[0]):
				crx = np.zeros((cells.shape[1]))
				cry = np.zeros((cells.shape[1]))
				for iv in range(cells.shape[1]):
					crx[iv] = nodes[int(cells[ic,iv]),0]
					cry[iv] = nodes[int(cells[ic,iv]),1]
				ax.plot(crx, cry, 'm-', linewidth = 0.5)

		else:
			raise ValueError("Mesh type not recognised!")

		ax.set_xlabel('$R \\; [m]$')
		ax.set_ylabel('$Z \\; [m]$')
		ax.set_title(this + " mesh")
		if top is not None: ax.set_ylim(top = top)
		if bottom is not None: ax.set_ylim(bottom = bottom)
		if left is not None: ax.set_xlim(left = left)
		if right is not None: ax.set_xlim(right = right)
		ax.set_aspect('equal')
	plt.show()

	return

#########################################################################################################

def plot_sv(where = ".", vmin = None, vmax = None, scale = "lin",
	        legend = False, loc = "best", grid = True, save = False):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	try:
		rates = pickle.load(open(os.path.join(where, 'rates.pkl'), 'rb'))
	except:
		rates = compute_rates(where = where)

	for key in rates.keys():

		if rates[key]['lnsv'] != []:

			(p0, p1, sv) = rates[key]['lnsv']
			sv = np.exp(sv) * 1E-06

			fig, ax = plt.subplots()

			for i1 in range(p1.shape[0]):
				try:
					if p1.min() < 1E+00:
						parameter = '$E_0 =$ '
						units     = ' $eV$' 
					else:
						parameter = '$n =$ '
						units     = ' $m^{-3}$' 
					if scale == 'lin': ax.semilogx(p0, sv[:,i1], label = parameter + '%.1E' % Decimal(str(p1[i1])) + units)
					else: ax.loglog(p0, sv[:,i1], label = parameter + '%.1E' % Decimal(str(p1[i1])) + units)
					legend = True
				except:
					if scale == 'lin': ax.semilogx(p0, sv)
					else: ax.loglog(p0, sv)
					legend = False
					break

			this  = ' '.join([rates[key]['database'], rates[key]['group'], rates[key]['reaction']])
			fancy_title = this + ' ' + '(' + rates[key]['type'] + ')' '\n' + rates[key]['equation']
			
			plt.xlabel(r'$T \; [eV]$')
			plt.ylabel(r"$\langle \sigma v \rangle \; [m^3 \cdot s^{-1}]$")
			if vmin is not None: plt.ylim(bottom = vmin)
			if vmax is not None: plt.ylim(top = vmax)
			plt.title(fancy_title)
			if legend is True: ax.legend(loc = loc)
			plt.grid(grid)
		
	plt.show()

	return


#########################################################################################################

def plot_rates(where = ".", database = None, group = None, reaction = None,
	           vmin = None, vmax = None, cmap = None, scale = "log", mask = True,
	   		   top = None, bottom = None, right = None, left = None, save = False):

	try:
		rates = pickle.load(open(os.path.join(where, 'rates.pkl'), 'rb'))
	except:
		rates = compute_rates(where = where)

	b2fgmtry = read_b2fgmtry(where = where, verbose = False, save = True)

	if isinstance(where, list) is True:
		raise ValueError('No multiple \"where\" accepted!')

	if not isinstance(vmin, list):  vmin  = [vmin]
	if not isinstance(vmax, list):  vmax  = [vmax]

	if database is None:
		database = []
		group    = []
		reaction = []
		for key in rates.keys():
			database.append(rates[key]['database'])
			group.append(rates[key]['group'])
			reaction.append(rates[key]['reaction'])
	else:
		if isinstance(database, list) is False: database = [database]
		if isinstance(group, list)    is False: group    = [group]
		if isinstance(reaction, list) is False: reaction = [reaction]

	for i in range(len(database)):

		try:    low = vmin[i]
		except: low = vmin[0]
		try:    up  = vmax[i]
		except: up  = vmax[0]

		for key in rates.keys():
			if (
				rates[key]['database'] == database[i] and
				rates[key]['group']    == group[i] and
				rates[key]['reaction'] == reaction[i]
				):
				try:
					value = rates[key]['reaction_rate']
					this  = ' '.join([database[i], group[i], reaction[i]])
					break
				# no data if database == ADAS or others
				except: pass

		if value != [] and value.max() > 0:

			if database[i] == 'AMMONX':
				triangles = read_triangle_mesh(where = where, verbose = False, save = False)
				cells = triangles['cells']
				nodes = triangles['nodes']
			else:
				(value, nodes, cells) = to_triangles(nsp = 1, vs = value, b2fgmtry = b2fgmtry)

			fancy_title = this + ' ' + '(' + rates[key]['type'] + ')' '\n' + rates[key]['equation']

			triplot(where = where, what = fancy_title, value = np.abs(value), sp = [''], nodes = nodes, cells = cells,
				    vmin = low, vmax = up, cmap = cmap, scale = scale, top = top, bottom = bottom,
				    right = right, left = left, mask = mask, save = save)
		else:
			print()
			print(this + ' missing...')

	plt.show()

	return

#########################################################################################################

def plot_mfps(where = ".", database = None, group = None, reaction = None,
	          vmin = None, vmax = None, cmap = None, scale = "log", mask = True,
	   		  top = None, bottom = None, right = None, left = None, save = False):

	try:
		rates = pickle.load(open(os.path.join(where, 'rates.pkl'), 'rb'))
	except:
		rates = compute_rates(where = where)

	b2fgmtry = read_b2fgmtry(where = where, verbose = False, save = True)

	if isinstance(where, list) is True:
		raise ValueError('No multiple \"where\" accepted!')

	if not isinstance(vmin, list):  vmin  = [vmin]
	if not isinstance(vmax, list):  vmax  = [vmax]

	if database is None:
		database = []
		group    = []
		reaction = []
		for key in rates.keys():
			database.append(rates[key]['database'])
			group.append(rates[key]['group'])
			reaction.append(rates[key]['reaction'])
	else:
		if isinstance(database, list) is False: database = [database]
		if isinstance(group, list)    is False: group    = [group]
		if isinstance(reaction, list) is False: reaction = [reaction]

	for i in range(len(database)):

		try:    low = vmin[i]
		except: low = vmin[0]
		try:    up  = vmax[i]
		except: up  = vmax[0]

		for mfp in ['0', '1']:

			for key in rates.keys():
				if (
					rates[key]['database'] == database[i] and
					rates[key]['group']    == group[i] and
					rates[key]['reaction'] == reaction[i]
					):
					try:
						value = rates[key]['mfp_' + mfp]
						this  = ' '.join([database[i], group[i], reaction[i]])
						break
					# no data if database == ADAS or others
					except: pass

			if value != [] and value.max() > 0:

				if database[i] == 'AMMONX':
					triangles = read_triangle_mesh(where = where, verbose = False, save = False)
					cells = triangles['cells']
					nodes = triangles['nodes']
				else:
					(value, nodes, cells) = to_triangles(nsp = 1, vs = value, b2fgmtry = b2fgmtry)

				fancy_title  = this + ' ' + '(' + rates[key]['type'] + ')' + ' [$' + rates[key]['species_' + mfp] + '$' + ' mfp]' + '\n' + rates[key]['equation']

				triplot(where = where, what = fancy_title, value = np.abs(value), sp = [''], nodes = nodes, cells = cells,
					    vmin = low, vmax = up, cmap = cmap, scale = scale, top = top, bottom = bottom,
					    right = right, left = left, mask = mask, save = save)
			else:
				if mfp == 0:
					print()
					print(this + ' missing...')

	plt.show()

	return
