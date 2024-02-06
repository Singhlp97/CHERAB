import numpy as np
import matplotlib.pyplot as plt
import math

from raysect.core import Point2D, Point3D, translate, Vector3D, rotate_basis
from raysect.optical import World, Spectrum
from raysect.primitive import Cylinder, Sphere
from raysect.optical.observer import PowerPipeline0D
from raysect.optical.observer.nonimaging.pixel import Pixel
from raysect.optical.material import AbsorbingSurface

from axisymmetric_mesh import axisymmetric_mesh_from_polygon
from axisymmetric_mesh_CUT_pieces import axisymmetric_mesh_from_polygon_CUT_pieces

from cherab.core.math.interpolators import Interpolate1DLinear, Interpolate2DLinear, Interpolate2DCubic

from my_functions.resize_contour_customized import resize_contour_customized
from my_functions.create_detectors import create_detectors
from my_functions.make_it_clockwise import make_it_clockwise

# External imports
import os
import time
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, pi

# Raysect imports
from raysect.core.workflow import SerialEngine
from raysect.core import translate, rotate_basis
from raysect.optical import World
from raysect.optical.observer.nonimaging.pixel import Pixel
from raysect.optical.observer import PowerPipeline0D
from raysect.primitive.mesh import Mesh
from raysect.optical.material.absorber import AbsorbingSurface

from raysect.primitive.mesh.obj import import_obj, export_obj
from raysect.primitive.mesh.vtk import export_vtk

########################################################################################################################
########################################################################################################################
########################################################################################################################

kind = "rear"

#in_directory = "radiation_load/input/Dpuff_Npuff_PFR.run.1.BGK=OFF.BCCON=8.2.5E20.Dpuff=1E21.Npuff=6E20.ANcoeff/SOLPS_data/"
#wall_points_files = ['baffle.extra']
# in_directory = "radiation_load/input/Dpuff_Npuff_PFR.run.1.BGK=OFF.BCCON=8.2.5E20.Dpuff=1E21.Npuff=6E20.ANcoeff.LJ-VA/mesh_wall_extra/"
# wall_points_files = ["hfs_plate_refined.extra"]
# in_directory = "/home/matteo.moscheni/raysect-cherab/runs/radiation_load/input/RF.mirror.VS_N_puff=6.00E19_Pe=2.27MW_Pi=1.49MW_sequence.run.13/mesh_observer_extra/"
# wall_points_files = ["1_HFStarg.txt"]
in_directory = "/home/lovepreet/CHERAB_files/radiation_261022/radiation_load/input/Fabio/mesh_absorbing_extra/"
wall_points_files = ["structure_c_.extra"]

out_directory = in_directory + "/../mesh_absorbing_vtk/"
os.system("mkdir " + out_directory)

index = 0

for wall_points_file in wall_points_files:

	with open(in_directory + wall_points_file, 'r') as fh:
		data = np.loadtxt(fh)

	basename = os.path.splitext(wall_points_file)[0]

	print()
	print(basename)
	print()

	data = data[::-1,:]

	wall_points_Nx2 = data

	#wall_points = wall_points_Nx2.copy()
	#wall_points = np.transpose(wall_points)

	#wall_points = make_it_clockwise(wall_points)

	wall_points_Nx2 = resize_contour_customized(np.transpose(wall_points_Nx2), 0.0)

	# wall_points_for_mesh must be in Nx2 format but resize_contour_customized returns 2xN 
	wall_points_for_mesh = wall_points_Nx2.T

	# wall_points_for_mesh = np.concatenate((wall_points_Nx2.T, np.zeros((1,2))))
	# wall_points_for_mesh = wall_points_for_mesh[::-1,:].T
	wall_points_for_mesh = np.array(wall_points_for_mesh)

	#wall_mesh = axisymmetric_mesh_from_polygon(wall_points_for_mesh)
	#export_vtk(wall_mesh, out_directory + './DEMO_wall_complete.vtk')

	wall_mesh = axisymmetric_mesh_from_polygon_CUT_pieces(wall_points_for_mesh)
	# wall_mesh = axisymmetric_mesh_from_polygon(wall_points_for_mesh)

	# file = open(out_directory + basename + '_refined.extra', 'w')
	# for i in range(wall_points_for_mesh.T.shape[0]):
	# 	file.write('{} {}\n'.format(wall_points_for_mesh.T[i,0], wall_points_for_mesh.T[i,1]))
	# file.close()
	export_vtk(wall_mesh, out_directory + basename + '_refined.vtk')
	export_obj(wall_mesh, out_directory + basename + '_refined.obj')

	index += 1
