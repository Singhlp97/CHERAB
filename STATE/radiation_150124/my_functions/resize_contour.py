""" author: Matteo Moscheni """

# resize of a given contour: the points of the latter are to be stored in a 2xN matrix
# if offset < 0 => reduced size
# if offset > 0 => enlarged size 

import numpy as np
import matplotlib.pyplot as plt

from raysect.core import Point2D, Point3D, Vector2D, Vector3D, rotate_basis

from my_functions.from_2xN_array_to_Point2D import from_2xN_array_to_Point2D
from my_functions.make_it_clockwise import make_it_clockwise


############
# function #

def resize_contour(contour_points_input, offset):

	contour_points = contour_points_input.copy()

	# to match with the definition of the function but not to need a change in the script
	# DO NOT EXAGGERATE either with positive and negative values: if too positive, you can get a contour 
	# that does no more fit with the overall geometry; if negative, strange switching between the position
	# of the translated points can appear.

	# error message
	if np.size(contour_points, 0) != 2: 
		print('\n\n***ERROR***\n\nThe points are to be stored in a narray with 2 rows and N columns!\n\n')

	SIZE_CONTOUR_POINTS = int(np.size(contour_points, 1))

	contour_points = make_it_clockwise(contour_points)

	#################################################################################################################
	#################################################################################################################

	#############################################
	# UNIFORMLY increasing the number of points #

	# the increase in the accuracy is meant as a decrease of the average distance between two consecutive points
	# in the mesh; it is obtained by splitting each reference segment in a certain number (factor_accuracy_UP - 1) of
	# sub-segments; in this way the initial number of points is enhanced by factor_accuracy_UP times; this procedure 
	# implies that the abovementioned distance is smaller for smaller segment and bigger for bigger ones but it
	# should be ok

	contour_points_new = []
	contour_points_dummy = []
	p_new = Point2D(0,0)

	accuracy_UP = input('\n\n 	Increase the number of points? (Y/N): ')

	if accuracy_UP == 'y' or accuracy_UP == 'Y':

	    factor_accuracy_UP = input('\n\n 	Enter the INTEGER increase factor (>1): ')
	    factor_accuracy_UP = int(factor_accuracy_UP)

	    for ii in range(0, SIZE_CONTOUR_POINTS):

	        p1x = contour_points[0][ii]
	        p1y = contour_points[1][ii]
	        p1 = Point2D(p1x, p1y)

	        # storage of the first point of the ii-th segment
	        contour_points_dummy = contour_points_dummy + [(p1)]

	        if ii < (SIZE_CONTOUR_POINTS - 1):
	            p2x = contour_points[0][ii + 1]
	            p2y = contour_points[1][ii + 1]
	        else :
	            p2x = contour_points[0][0]
	            p2y = contour_points[1][0]

	        p2 = Point2D(p2x, p2y)   

	        y_vector = p1.vector_to(p2)

	        # new points creation: the vector between the extrema p1 and p2 is divided in factor_accuracy_UP parts
	        # so that to identify the (factor_accuracy_UP - 1) inner points
	        y_vector = y_vector / factor_accuracy_UP

	        for jj in range(1, factor_accuracy_UP):

	            p_new = p1 + y_vector * jj
	            # storage of the jj-th new point of the ii-th segment
	            contour_points_dummy = contour_points_dummy + [(p_new)]

	        # contour_points_dummy = contour_points_dummy + [(p2)]
	        # NO NEED: it will be automatically inserted as the p1 at the step (ii + 1)-th


	    # the previous vector is overwritten
	    contour_points_new = contour_points_dummy.copy()
	    # to be updated! The new one is factor_accuracy_UP times the old one
	    SIZE_CONTOUR_POINTS = int(np.size(contour_points_new)) 

	elif accuracy_UP == 'n' or accuracy_UP == 'N':
		# contour_points_new needs to be a tuple with Point2D as elements
		# => trnasforming the 2xN-array contour_points in such an object while leaving it unchanged
		# since no accuracy_UP is required
		contour_points_new = from_2xN_array_to_Point2D(contour_points)

	else :
	    print('\n\n***ERROR***\n\nInvalid input! The contour points are taken as reference with no increase!')
	    contour_points_new = from_2xN_array_to_Point2D(contour_points)


	#################################################################################################################
	#################################################################################################################

	##########
	# resize #

	# the set of mesh_points is traveled CLOCKWISELY
	# it is chosen to have OUTWARD-POINTING normal unit vectors for the mesh points

	contour_points_resized = []

	y_vector_1 = []
	y_vector_2 = []
	y_normal_1 = []
	y_normal_2 = []

	y_vector = []
	y_normal = []

	Y_NORMAL = []
	NEW_POINT = []

	# contour_points_new is now an array of Point2D!!!

	# loop to create the new set of inner points of the detectors
	for ii in range(0, SIZE_CONTOUR_POINTS):

	    # if ii = 0 -> contour_points(-1) = the last one OK
	    p1x = contour_points_new[ii - 1].x
	    p1y = contour_points_new[ii - 1].y
	    p1 = Point2D(p1x, p1y)

	    p2x = contour_points_new[ii].x
	    p2y = contour_points_new[ii].y
	    p2 = Point2D(p2x, p2y)

	    if ii < (SIZE_CONTOUR_POINTS - 1): 
	        p3x = contour_points_new[ii + 1].x
	        p3y = contour_points_new[ii + 1].y
	    else :
	        p3x = contour_points_new[0].x
	        p3y = contour_points_new[0].y

	    p3 = Point2D(p3x, p3y)

	    # take the normal between p1 and p2 and the one between p2 and p3: they have the same magnitude and hence 
	    # the rule of the sum of two vectors gives back the vector that bisects the angle formed by thw two of them;
	    # consequently, switching the sign, for sure the resulting vector points in a direction that does not
	    # intersect any other segment; at a distance of DETECTOR_OFFSET along this direction, the new inner point
	    # (edge of the detector) is placed

	    # evaluate y_vector
	    y_vector_1 = p1.vector_to(p2)
	    y_vector_2 = p2.vector_to(p3)
	        
	    # normal vectors: rotation of the y_vectors in such a way that the normal ones are outward-pointing
	    y_normal_1 = Vector2D(- y_vector_1.y, y_vector_1.x).normalise()
	    y_normal_2 = Vector2D(- y_vector_2.y, y_vector_2.x).normalise()

	    # evaluate the "average normal"
	    Y_NORMAL = (y_normal_1 + y_normal_2).normalise()
	    NEW_POINT = p2 + Y_NORMAL * offset

	    contour_points_resized = contour_points_resized + [(NEW_POINT.x, NEW_POINT.y)]

	# plt.figure()
	# for ii in range(np.size(contour_points,1)-1):
	# 	plt.plot([contour_points[0,ii], contour_points[0,ii+1]], [contour_points[1,ii],contour_points[1,ii+1]], 'k-')
	# 	plt.plot([contour_points_resized[ii][0], contour_points_resized[ii+1][0]], [contour_points_resized[ii][1],contour_points_resized[ii+1][1]], 'm-')
	# plt.axis('equal')
	# plt.show()

	contour_points_resized_ok = np.transpose(contour_points_resized)

	return contour_points_resized_ok
