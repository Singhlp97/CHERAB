import numpy as np
import matplotlib.pyplot as plt
import math

from raysect.core import Point2D, Point3D, translate, Vector3D, rotate_basis

def from_Point2D_to_2xN_array(points_Point2D_input):

	points_Point2D = points_Point2D_input.copy()

	NUM_POINTS = np.size(points_Point2D)
	points_2xN_array = np.zeros((2, NUM_POINTS))

	for ii in range(0, NUM_POINTS):
		points_2xN_array[0][ii] = points_Point2D[ii].x
		points_2xN_array[1][ii] = points_Point2D[ii].y

	return points_2xN_array

