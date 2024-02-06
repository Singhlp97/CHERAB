""" author: Matteo Moscheni """

# generate the detectors starting from the 2xN narray of the coordinates of their edges

import numpy as np
import matplotlib.pyplot as plt

from raysect.core import Point2D, Point3D, translate, Vector2D, Vector3D, rotate_basis
from raysect.optical import World, Spectrum
from raysect.primitive import Cylinder
from raysect.optical.observer import PowerPipeline0D
from raysect.optical.observer.nonimaging.pixel import Pixel
from raysect.optical.material import AbsorbingSurface

from cherab.core.math import sample2d, AxisymmetricMapper
#from cherab.tools.emitters import RadiationFunction
#from cherab.tools.primitives import axisymmetric_mesh_from_polygon
from my_functions.make_it_clockwise import make_it_clockwise

def create_detectors(detector_points_input):

    detector_points = detector_points_input.copy()
    detector_points = make_it_clockwise(detector_points)

    x_width = 1e-2

    p1x = 0
    p1y = 0
    p1 = []

    p2x = 0
    p2y = 0
    p2 = []

    y_vector = []
    y_vector_full = []
    y_normal = []
    y_width = 0
    detector_center = []
    detectors = []

    SIZE_DETECTORS_POINTS = np.size(detector_points, 1)

    # now they are known the edge points of each detector and it is then possible to construct them in the usual way
    # notice that it is NOT possible to add this procedure in the previous loop because it is here needed the 
    # simultaneus knowledge of both the edges of the detectors but here above they are computed one by one
    for ii in range(0, SIZE_DETECTORS_POINTS):

        p1x = detector_points[0][ii]
        p1y = detector_points[1][ii]
        p1 = Point3D(p1x, 0, p1y)

        # alternatively you can start with p1 = ... [-1] and no if/else will be needed
        if ii < (SIZE_DETECTORS_POINTS - 1):
            p2x = detector_points[0][ii + 1]
            p2y = detector_points[1][ii + 1]
        else :
            p2x = detector_points[0][0]
            p2y = detector_points[1][0]

        p2 = Point3D(p2x, 0, p2y)

        y_vector_full = p1.vector_to(p2)
        y_vector = y_vector_full.normalise()
        y_normal = Vector3D(y_vector.z, 0, - y_vector.x).normalise() # inward pointing

        # evaluate the central point of the detector
        detector_center = p1 + y_vector_full * 0.5
        # evaluate the length in the poloidal section
        y_width = y_vector_full.length

        # to populate it step by step: added p1 and p2 in the end
        detectors = detectors + [(ii, x_width, y_width, detector_center, y_normal, y_vector)]

        
        if ii < SIZE_DETECTORS_POINTS:
            plt.plot([p1x, p2x], [p1y, p2y], 'm')
            plt.plot([p1x, p2x], [p1y, p2y], '.m')
            pc = detector_center
            pcn = pc + y_normal * 0.1
            plt.plot([pc.x, pcn.x], [pc.z, pcn.z], 'k')
    
    plt.axis('equal')
        

    return detectors
