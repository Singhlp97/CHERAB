import numpy as np
import matplotlib.pyplot as plt

from raysect.core import Point2D
from cherab.core.math import PolygonMask2D

####################################################################################################################################################

# constants

me = 9.11E-31 # [kg]
mp = 1.67E-27 # [kg]

####################################################################################################################################################

def plot_1d(x = None, y = None, y0 = None, xlabel = None, ylabel = None, scale = "lin", title = None, output = False):

	# plt.rcParams.update({'font.size': 14})
	linewidth = 3

	if type(y) is tuple:
		num_curves = 2
		ytuple = y
		yraws = []
	else:
		num_curves = 1
		ytuple = (y,)

	for y in ytuple:

		if type(y) is not np.ndarray:
			(there_is_wall, wall, inside_wall) = make_chamber_wall()
			num = 500
			x = np.linspace(x.min(), x.max(), num)
			y_ = np.zeros((np.size(x)))
			for j in range(y_.shape[0]):
				y_[j] = y(x[j], y0)
				if there_is_wall:
					y_[j] *= inside_wall(x[j], y0)
		else:
			if num_curves == 1:
				y_ = y
			else:
				yraws += [y]

		if num_curves == 1:
			plt.figure()
			fig = plt.gcf()
			fig.set_size_inches(7, 5.5)
			if   scale == "lin": plt.plot(x, y_, linewidth = linewidth)
			elif scale == "log": plt.semilogy(x, y_, linewidth = linewidth)
			plt.xlabel(xlabel)
			plt.ylabel(ylabel)
			plt.title(title)

	if num_curves == 2:
		fig, ax1 = plt.subplots(figsize = (7, 5.5))
		ax2 = ax1.twinx()
		ax1.plot(x, yraws[0], 'b-', linewidth = linewidth)
		ax2.plot(x, yraws[1], 'm-', linewidth = linewidth)
		ax1.set_xlabel(xlabel)
		ax1.set_ylabel(ylabel[0], color = 'b')
		ax2.set_ylabel(ylabel[1], color = 'm')
		ax1.set_title(title)

	if output is True: return (x, y_)
	else: return

####################################################################################################################################################

def plot_2d(ri = None, zi = None, f2d = None, scale = None, vmin = None, vmax = None, eq = None, wall = None, title = None, cmap_label = None):

	(there_is_wall, wall, inside_wall) = make_chamber_wall()

	if type(f2d) is not np.ndarray:
		num = 500
		ri = np.linspace(ri.min(), ri.max(), num)
		zi = np.linspace(zi.min(), zi.max(), num)
		data = np.zeros((np.size(zi), np.size(ri)))
		for i in range(data.shape[0]):
			for j in range(data.shape[1]):
				data[i,j] = f2d(ri[j], zi[i])
				if there_is_wall:
					data[i,j] *= inside_wall(ri[j], zi[i])
	else:
		data = f2d

	data[np.where(data == 0)] = np.nan
	if vmin is not None: data[np.where(data < vmin)] = vmin
	if vmax is not None: data[np.where(data > vmax)] = vmax
	if scale == "log":   data = np.log10(data)

	fig, ax = plt.subplots()
	c = ax.pcolormesh(ri, zi, data, cmap = 'jet')
	fig.colorbar(c, ax = ax, label = cmap_label)
	if eq is not None:
		xlcfs = np.append(eq.lcfs_polygon[:,0], np.array(eq.lcfs_polygon[0,0]))
		ylcfs = np.append(eq.lcfs_polygon[:,1], np.array(eq.lcfs_polygon[0,1]))
		plt.plot(xlcfs, ylcfs, 'm-')
	if wall is not None:
		plt.plot(wall[:,0], wall[:,1], 'k-', linewidth = 1)
	ax.set_xlabel('$R \\; [m]$')
	ax.set_ylabel('$Z \\; [m]$')
	# ax.set(xlim = (0, 1), ylim = (-1, 1))
	# ax.set_aspect('equal', adjustable = 'box')
	ax.set_title(title)
	ax.axis('equal')

####################################################################################################################################################

def normalise_psi(eq = None):
		return (eq.psi_data - eq.psi_axis) / (eq.psi_lcfs - eq.psi_axis)

####################################################################################################################################################

def ExponentialDecay(x = None, x0 = None, decay_length = None):
	return np.exp(- (x - x0) / decay_length)

####################################################################################################################################################

def Gaussian(x = None, x0 = None, smoothing = 1.0):
	return np.exp(-(x - x0)**2 / smoothing)

####################################################################################################################################################

def FermiDirac(psi = None, psi_lcfs = None, smoothing = 1.0):
	return 1 / (1 + np.exp((psi_lcfs - psi) / smoothing))

####################################################################################################################################################

def SoundSpeed(Ti = None, Te = None, mi = 2 * mp, adiabatic_coefficient = 1.0):
	return np.sqrt((adiabatic_coefficient * Ti + Te) / mi)

####################################################################################################################################################

def DecayLengthStangeby(Ti = None, Te = None, mi = 2 * mp, diffusivity = None, connection_length = None):
	return np.sqrt(diffusivity * connection_length / SoundSpeed(Ti = Ti, Te = Te, mi = mi))

####################################################################################################################################################

def DensityDoublingStangeby(x = None, y = None, eq = None):
	if eq.inside_lcfs(x,y) == 0:
		centre = eq.magnetic_axis
		point  = Point2D(x,y)
		vector = centre.vector_to(point).normalise()
		return (1.0 + np.arccos(vector.x) / np.pi)
	else: return 1.0

####################################################################################################################################################

def ComputeConnectionLength(eq = None):

	# - interpolate B map from EFIT
	# - build axisymmetric 3D equilibrium
	# - integrate along field line ??? psi = constant?

	return 

####################################################################################################################################################

def ComputeDecayLength(R = None, data = None, datamin = None, title = None):

	Rlcfs = 0.6916 # [m] at OMP Z = 0

	ok = data > datamin
	sol = R > Rlcfs
	sol *= ok

	Rsol = R[sol]
	datasol = data[sol]

	(x1, y1) = (Rsol[0], np.log(datasol[0]))
	(x2, y2) = (Rsol[-1], np.log(datasol[-1]))

	ldecay = np.abs((x2 - x1) / (y2 - y1))

	plt.figure()
	plt.semilogy(R, data, '-', linewidth = 4)
	plt.semilogy(Rsol, datasol[0] * np.exp( - (Rsol - Rlcfs) / ldecay), '--', linewidth = 3)
	plt.title(title + ' $\\lambda = $ {:.4G} m'.format(ldecay))

####################################################################################################################################################

def make_chamber_wall(verbose = False):

	try:
		wall = np.loadtxt('./mesh.extra')
		inside_wall = PolygonMask2D(wall[:-1,:2])
		there_is_wall = True
		if verbose is True: print('Wall loaded from ./mesh.extra...')
	except:
		wall = None
		inside_wall = None
		there_is_wall = False
		if verbose is True: print('No wall found...')

	return (there_is_wall, wall, inside_wall)