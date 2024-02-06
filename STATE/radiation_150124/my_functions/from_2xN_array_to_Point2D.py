import numpy as np
import math

from raysect.core import Point2D

def from_2xN_array_to_Point2D(array_2xN_input):

	array_2xN = array_2xN_input.copy()

	NN = np.size(array_2xN, 1)

	points_Point2D = []

	for ii in range(0, NN):
		pp = Point2D(0,0)
		pp.x = array_2xN[0][ii]
		pp.y = array_2xN[1][ii]
		points_Point2D = points_Point2D + [(pp)]

	return points_Point2D



