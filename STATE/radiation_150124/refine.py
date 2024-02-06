import os
import numpy as np
from raysect.primitive.mesh.obj import export_obj
from raysect.primitive.mesh.vtk import export_vtk
from raysect.core import Point2D, Vector2D

from axisymmetric_mesh import axisymmetric_mesh_from_polygon
from axisymmetric_mesh_CUT_pieces import axisymmetric_mesh_from_polygon_CUT_pieces


in_directory = "radiation_load/input/Dpuff_Npuff_PFR.run.1.BGK=OFF.BCCON=8.2.5E20.Dpuff=1E21.Npuff=6E20.ANcoeff/SOLPS_data/"
wall_points_files = ['mesh.extra']

out_directory = in_directory + "/../mesh_wall_vtk/"
os.system("mkdir " + out_directory)

index = 0

for wall_points_file in wall_points_files:

	with open(in_directory + wall_points_file, 'r') as fh:
		data = np.loadtxt(fh)

	basename = os.path.splitext(wall_points_file)[0]

	print()
	print(basename)
	print()

	#data = data[::-1,:]

	wall_points_Nx2 = np.zeros((np.size(data,0), 2))

	wall_points_Nx2[:,0] = data[:,0]
	wall_points_Nx2[:,1] = data[:,1]

	pts = []

	for (i,pt) in enumerate(wall_points_Nx2):

		if i < wall_points_Nx2.shape[0]-1:
			num = 1
			p1 = Point2D(pt[0], pt[1])
			p2 = Point2D(wall_points_Nx2[i+1,0], wall_points_Nx2[i+1,1])
			l12 = p1.distance_to(p2)
			dl = l12 / num
			v12 = p1.vector_to(p2).normalise()

			for i in range(num-1):
				p = p1 + v12 * dl * i
				pts.append([p.x, p.y])

	pts = np.array(pts)
	#wall_mesh = axisymmetric_mesh_from_polygon_CUT_pieces(pts)
	wall_mesh = axisymmetric_mesh_from_polygon(data)
	export_vtk(wall_mesh, out_directory + basename + '.vtk')
	export_obj(wall_mesh, out_directory + basename + '.obj')

	index += 1
