import os
import json
import numpy as np
import matplotlib.pyplot as plt

from raysect.core import Point3D, Vector3D

def CreateSmartPixels(cfg = None, observation_line = None, observation_point = None):

	"""

	TODO MMM:
	- save figure of pixel grid in 2D plane
	
	"""

	###############

	# CAUTION.
	#
	# reference simulation must be of same "type" w.r.t. current:
	# - same input data
	# - same type_radiation
	# - same step, astra, b2, eirene, etc. switched (True/False)
	# ...
	#
	# under user-control:
	# - specific run number 
	# - .csv filename

	# current cfg

	basedir = os.path.join(cfg['baserun'], 'output', cfg['run'], cfg['plasma']['type_radiation'])
	case    = cfg['raytracing']['observer']['specs']['smart_pixelling']['case']

	plot_smart_pixelling      = cfg['plotting']['plot_smart_pixelling']
	save_figures              = cfg['plotting']['save_figures']
	output_directory_extended = cfg['output_directory_extended']

	del cfg

	# cfg of reference simulation

	# load .json stored in output directory
	cfg = json.load(open(os.path.join(basedir, case, 'input', 'configFile.json'), 'r'))
	# load data in .csv
	data = np.loadtxt(open(os.path.join(basedir, case, 'camera_shift.csv'), 'r'), delimiter = ',', skiprows = 0)
	###############

	# sum through columns and rows

	# watch out index convention
	fsumy = data.sum(axis = 0)
	fsumx = data.sum(axis = 1)

	Nx = cfg['raytracing']['observer']['specs']['nx_pixels']
	Ny = cfg['raytracing']['observer']['specs']['ny_pixels']

	nx = Nx * fsumx / fsumx.sum()
	ny = Ny * fsumy / fsumy.sum()

	###############

	# start building uniformly-distributed pixel grid (i.e. original ones)

	# CAUTION.
	#
	# PinholeCamera: rays are fired from observation_point
	#
	# VectorCamera: rays are fired from pixel centres!
	# Location of image plane hence becomes important!

	fov = cfg['raytracing']['observer']['specs']['fov']
	offset = cfg['raytracing']['observer']['surface_offset']
	image_width = 2 * np.tan(np.pi / 180 * 0.5 * fov) * offset

	#observation_point -= observation_line.normalise() * 0.5
	observation_point = Point3D(observation_point.x,
								observation_point.y,
								observation_point.z)

	# image centre offset-m away from observation_point
	image_centre = observation_point + offset * observation_line
	y_direction = Vector3D(0, 0, 1)
	x_direction = y_direction.cross(observation_line).normalise()

	y_start = image_centre - (x_direction + y_direction) * 0.5 * image_width
	y_end = y_start + y_direction * image_width

	x_start = image_centre - (x_direction + y_direction) * 0.5 * image_width
	x_end = x_start + x_direction * image_width

	###############

	# build pixel centres along x

	dlx = image_width / Nx
	px = [Point3D(x_start.x, x_start.y, x_start.z)]

	px_smart = []
	nxs = []

	for i in range(Nx-1):

		p = x_start + (i+1) * dlx * x_direction

		px.append(Point3D(p.x, p.y, p.z))

		n = np.max([int(np.ceil(np.max([nx[i], nx[i+1]]))), 1])

		nxs.append(n)

		dl = px[i].distance_to(px[i+1]) / n

		for j in range(n):
			p = px[i] + j * dl * x_direction
			px_smart.append(Point3D(p.x, p.y, p.z))

	###############

	# build pixel centres along y

	dly = image_width / Ny
	py = [Point3D(y_start.x, y_start.y, y_start.z)]

	py_smart = []
	nys = []

	for i in range(Ny-1):

		p = y_start + (i+1) * dly * y_direction

		py.append(Point3D(p.x, p.y, p.z))

		n = np.max([int(np.ceil(np.max([ny[i], ny[i+1]]))), 1])

		nys.append(n)

		dl = py[i].distance_to(py[i+1]) / n

		for j in range(n):
			p = py[i] + j * dl * y_direction
			py_smart.append(Point3D(p.x, p.y, p.z))

	###############

	# build pixel_origins and pixel_directions
	# ~ lexicographic representation

	pixel_origins    = []
	pixel_directions = []

	plot_x = []
	plot_y = []
	plot_z = []

	#plt.figure()
	#ax = plt.axes(projection='3d')

	for ix in range(len(px_smart)):

		tmp_origins    = []
		tmp_directions = []

		for iy in range(len(py_smart)):

			ly = py[0].distance_to(py_smart[iy])
			p = px_smart[ix] + ly * y_direction
			p = Point3D(p.x, p.y, p.z)
			tmp_origins.append(p)
			tmp_directions.append(observation_point.vector_to(p))

			#ax.plot3D([observation_point.x, p.x],
			#		  [observation_point.y, p.y],
			#		  [observation_point.z, p.z])

			plot_x.append(p.x)
			plot_y.append(p.y)
			plot_z.append(p.z)

		pixel_origins    += [tmp_origins]
		pixel_directions += [tmp_directions]

	###############

	print()
	print('Number original pixels along x: {}'.format(Nx))
	print('Number smart    pixels along x: {}'.format(len(px_smart)))
	print()
	print('Number original pixels along y: {}'.format(Ny))
	print('Number smart    pixels along y: {}'.format(len(py_smart)))
	print()

	if plot_smart_pixelling is True:

		#plt.figure()
		#ax = plt.axes(projection='3d')
		#ax.scatter3D(plot_x, plot_y, plot_z)

		plt.figure()
		plt.plot(np.arange(0, Ny), ny)
		plt.plot(np.arange(0, Ny-1), nys, 's-')
		plt.xlabel('pixel number')
		plt.ylabel('number of smart pixels per default pixel')
		plt.title('horizontal pixels')

		if save_figures is True: plt.savefig(os.path.join(output_directory_extended, 'py.png'))

		plt.figure()
		plt.plot(nx, np.arange(0, Nx))
		plt.plot(nxs, np.arange(0, Nx-1), 's-')
		plt.xlabel('number of smart pixels per default pixel')
		plt.ylabel('pixel number')
		plt.title('vertical pixels')

		if save_figures is True: plt.savefig(os.path.join(output_directory_extended, 'px.png'))

		plt.show()

	return (pixel_origins, pixel_directions)