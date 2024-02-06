import os
import numpy as np
import pickle

from solps_python_scripts.read_ft33 import read_ft33
from solps_python_scripts.read_ft34 import read_ft34
from solps_python_scripts.read_ft35 import read_ft35

def read_triangle_mesh(where = ".", verbose = True, save = None):

	# triangles = read_triangle_mesh(fort33,fort34,fort35)
	#
	# Wrapper routine to read all triangle data at once.
	#
	# Returns nodes, cells, nghbr, side and cont as fiels of triangles-struct.


	# Author: Wouter Dekeyser
	# November 2016
	#
	# Re-writte in python by: Matteo Moscheni
	# E-mail: matteo.moscheni@tokamakenergy.co.uk
	# February 2022

	"""
		TODO mmm:
		- include fort.3* in read_ft3*.py already
	"""

	triangles = {}

	triangles['nodes'] = read_ft33(where = where, verbose = verbose)
	triangles['cells'] = read_ft34(where = where)
	links              = read_ft35(where = where)

	triangles['nghbr'] = links['nghbr']
	triangles['side']  = links['side']
	triangles['cont']  = links['cont']
	triangles['ixiy']  = links['ixiy']

	triangles['mesh_extent'] = {}

	triangles['mesh_extent']['minr'] = triangles['nodes'][:,0].min()
	triangles['mesh_extent']['maxr'] = triangles['nodes'][:,0].max()
	triangles['mesh_extent']['minz'] = triangles['nodes'][:,1].min()
	triangles['mesh_extent']['maxz'] = triangles['nodes'][:,1].max()

	if save is True:
		pickle.dump(triangles, open(os.path.join(where, "triangle_mesh.pkl"), "wb"))

	return triangles
