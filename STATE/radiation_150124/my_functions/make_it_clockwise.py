""" author: Matteo Moscheni """

# the routine takes a closed contour and, if not, make it be traveled clockwisely.
# the contour is to be a 2xN array

import numpy as np
import matplotlib.pyplot as plt

from raysect.core import Point2D, Point3D, Vector2D, Vector3D, rotate_basis

from my_functions.from_2xN_array_to_Point2D import from_2xN_array_to_Point2D

############
# function #

def make_it_clockwise(contour_points_input):

	contour_points = contour_points_input.copy()
	SIZE_CONTOUR_POINTS = np.size(contour_points, 1)

	if SIZE_CONTOUR_POINTS == 1:
		print('\n\n***ERROR***\n\nThe contour is a single point! Not acceptable...')
		return 0

	#####################
	# check orientation #

	# ANALYTICAL CHECK that points are orientated such that to travel the closed path CLOCKWISELY

	# initializations
	point_REF = Point2D(0,0)
	vec = Vector2D(0,0)
	sum_vec = Vector2D(0,0)

	# evaluation of the internal, reference node as the average of the points (VECTORS) of the mesh:
	# vecAvg = (vec1 + vec2 + vec3 + ...) / number_of_vectors (center-of-mass like procedure with
	# a fixed weigth <-> mass)

	# ***ATTENTION***
	# there are cases (e.g. very thin but distorted shapes) in which this does not work but with
	# delta_2pi it is reduced the probability of not noticing this evenience.
	for ii in range(0, SIZE_CONTOUR_POINTS):
	    vec.x = contour_points[0][ii]
	    vec.y = contour_points[1][ii]
	    sum_vec = sum_vec + vec

	vec = sum_vec / SIZE_CONTOUR_POINTS
	# reference inner point
	point_REF.x = vec.x
	point_REF.y = vec.y

	# initializations
	v1x = 0
	v1y = 0
	v2x = 0
	v2y = 0
	v1 = Vector2D(0,0)
	v2 = Vector2D(0,0)
	v_cross_z = 0
	FLAG = 0

	num_loops = 0

	# FLAG = 1 -> path traveled clockwisely
	# FLAG = 0 -> path traveled counterclockwisely -> NO
	# while loop => if the iterations are more than one it means that something is wrong with the contour
	while FLAG == 0:

		num_loops = num_loops + 1
		sum_angle = 0

		for ii in range(0, SIZE_CONTOUR_POINTS):
			v1x = contour_points[0][ii] - point_REF.x
			v1y = contour_points[1][ii] - point_REF.y
			v1 = Vector2D(v1x, v1y)

			if ii < (SIZE_CONTOUR_POINTS - 1):
				v2x = contour_points[0][ii + 1] - point_REF.x
				v2y = contour_points[1][ii + 1] - point_REF.y
			else :
				v2x = contour_points[0][0] - point_REF.x
				v2y = contour_points[1][0] - point_REF.y

			v2 = Vector2D(v2x, v2y)

			# ***ATTENTION***
	        # to evaluate the angle theta span while traveled along the contour, it is needed to translate 
	        # the contour in the baricentric point (point_REF). However, if such a point coincides with
	        # a point of the contour, it is met a point in zero that is not allowed in the expression
	        # abs(np.arccos(v1.dot(v2)/v1.length/v2.length)): 0/0 that is, in this case, 0.
	        # sum_angle is the sum of the angle-with-sign formed by one vector with the adjacent one;
	        # the sign is given by the cross product: positive if from vec1 to vec2 there is a counterclockwise
	        # rotation, negative viceversa -> the sum is to be (-1) * 2 * pi approximately;
	        # theta between vec1 and vec2 is evaluated back from the definition of the dot product

			if v1.length != 0 and v2.length != 0:
				v_cross_z = v1.cross(v2)
				sum_angle = sum_angle + np.sign(v_cross_z) * abs(np.arccos(v1.dot(v2)/v1.length/v2.length))
			# else: sum_angle = sum_angle + 0

		#if sum_angle >= (-2*np.pi - delta_2pi) and sum_angle <= (-2*np.pi + delta_2pi) :
		if num_loops <= 2 and sum_angle < 0:
			FLAG = 1
		elif num_loops <= 2 and sum_angle > 0:
			#elif sum_angle >= (2*np.pi - delta_2pi) and sum_angle <= (2*np.pi + delta_2pi) :
			print('\n\n***ATTENTION*** The path is NOT traveled clockwisely! The set of points is flipped...\n\n')
			# in this case the set of points is flipped and re-checked
#			contour_points[0][:] = contour_points[0][::-1]
#			contour_points[1][:] = contour_points[1][::-1]
		else:
			print('\n\n***ERROR***\n\nThere is something wrong with the contour! It is left unchanged...\n\n')
			break


	return contour_points
