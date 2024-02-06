""" author: Matteo Moscheni """

# source_points_input = 2xN array with the location in the poloidal section of the points of emission
# power_profile_1D_input = 1xN array with the power per unit volume emitted in each of the previous points

import numpy as np
import matplotlib.pyplot as plt
import math

from raysect.core import Point2D, Point3D, translate, Vector3D, rotate_basis
from raysect.optical import World, Spectrum
from raysect.primitive import Cylinder
from raysect.optical.observer import PowerPipeline0D
from raysect.optical.observer.nonimaging.pixel import Pixel
from raysect.optical.material import AbsorbingSurface

from cherab.core.math import sample2d, AxisymmetricMapper
from cherab.tools.emitters import RadiationFunction
from cherab.tools.primitives import axisymmetric_mesh_from_polygon

from cherab.core.math.interpolators import Interpolate2DLinear, Interpolate2DCubic

from my_functions.regular_polygon import regular_polygon
from my_functions.resize_contour import resize_contour
from my_functions.create_detectors import create_detectors
from my_functions.from_Point2D_to_2xN_array import from_Point2D_to_2xN_array
from my_functions.progress_bar_Ricky import progress_bar_Ricky

def make_them_ordered_my_way(source_points_input):

	source_points = source_points_input.copy()
	num_source_points = np.size(source_points, 1)

	source_points_ordered = np.zeros((2, num_source_points))

	index_min_ordinate = 0

	for ii in range(num_source_points):

		min_ordinate = math.inf

		for jj in range(np.size(source_points, 1)):
 
			if source_points[1][jj] < min_ordinate:
				min_ordinate = source_points[1][jj]
				index_min_ordinate = jj

			elif source_points[1][jj] == min_ordinate:
				if source_points[0][jj] < source_points[0][index_min_ordinate]:
					min_ordinate = source_points[1][jj]
					index_min_ordinate = jj

		source_points_ordered[0][ii] = source_points[0][index_min_ordinate]
		source_points_ordered[1][ii] = source_points[1][index_min_ordinate]

		source_points[1][index_min_ordinate] = math.inf

	return source_points_ordered

