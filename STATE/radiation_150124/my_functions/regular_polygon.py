""" author: Matteo Moscheni """

# regular polygon generator: given the central point (Point2D), the radius and the number of sides, 
# a regular polygon is returned

import numpy as np
import matplotlib.pyplot as plt
import math

from raysect.core import Point2D

from my_functions.make_it_clockwise import make_it_clockwise

############
# function #

def regular_polygon(centre_point, radius, number_sides):

	# internal angle to displace the points along the circumference
	d_angle = (2*np.pi) / number_sides

	pr = []
	pz = []
	polygon_points = np.zeros((2, number_sides))

	for i in range(number_sides):
		# from the central point it is evaluated the new one on the circumference
	    pr = centre_point.x + radius * np.sin(d_angle * i)
	    pz = centre_point.y + radius * np.cos(d_angle * i)
	    polygon_points[:, i] = pr, pz

	polygon_points = make_it_clockwise(polygon_points)

	return polygon_points

