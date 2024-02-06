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

from my_functions.from_2xN_array_to_Point2D import from_2xN_array_to_Point2D
from my_functions.from_Point2D_to_2xN_array import from_Point2D_to_2xN_array
from my_functions.make_it_clockwise import make_it_clockwise

def resize_contour_customized(contour_points_input, offset):

    mesh_points = contour_points_input.copy()
    mesh_points = make_it_clockwise(mesh_points)

    SIZE_MESH_POINTS = np.size(mesh_points, 1)

    # Plot of the reference contour with labels
    fig_1 = plt.figure(1)
    ax = fig_1.add_subplot(111)
    plt.plot(mesh_points[0,:], mesh_points[1,:], 'ko-')
    kk = 0
    for x,y in zip(mesh_points[0,:], mesh_points[1,:]):
        label = "{}".format(kk)
        kk = kk + 1
        plt.annotate(label, # this is the text
                     (x,y), # this is the point to label
                     textcoords="offset points", # how to position the text
                     xytext=(0,8), # distance from text to points (x,y)
                     ha='left') # horizontal alignment can be left, right or center

    plt.plot([mesh_points[0, -1], mesh_points[0, 0]], [mesh_points[1, -1], mesh_points[1, 0]], 'ko-')
    plt.title('Reference contour')
    plt.axis('equal')
    plt.show()

    mesh_points = from_2xN_array_to_Point2D(mesh_points)

    #################################################################################################################
    #################################################################################################################

    ###########################################
    # NON-UNIFORMLY increasing the resolution #

    # the increase in the accuracy is meant as a decrease of the average distance between two consecutive points
    # in the contour; it is obtained by splitting each reference segment in a certain number 
    # (factor_accuracy_UP - 1) of sub-segments; 
    # the procedure is implemented in a loop so that it is possible to customize the contour as much as possible,
    # e.g. enhancing of a certain factor in a certain range and of another factor in another portion.

    # initialization to enter the loop
    accuracy_UP = 'y'
    num_accuracy_UP = 0

    while accuracy_UP == 'y' or accuracy_UP == 'Y':

        if num_accuracy_UP == 0:
            accuracy_UP = input('\n\n######################################\n\n*** INCREASE THE RESOLUTION SOMEWHERE? (Y/N): ')
        else:
            accuracy_UP = input('\n\n######################################\n\n*** FURTHER INCREASE THE RESOLUTION SOMEWHERE? (Y/N): ')

        # YES -> the resolution is to be enhanced somewhere
        # NO -> the resolution is fine with the points of the contour
        if accuracy_UP == 'y' or accuracy_UP == 'Y':

            num_accuracy_UP = num_accuracy_UP + 1

            # acquisition of the range in whcih to increase the resolution
            print('\n\nIn which range do you need to increase the resolution? Notice that 0 => 0 represents the entire contour')
            where_accuracy_UP_1 = input('\nEnter the first index (counterclockwisely, min = 0 and max = {}): '.format(SIZE_MESH_POINTS - 1))
            where_accuracy_UP_2 = input('\nEnter the second index (counterclockwisely, min = 1 and max = 0): ')

            # they need to be integers since they are indeces
            where_accuracy_UP_1 = int(where_accuracy_UP_1)
            where_accuracy_UP_2 = int(where_accuracy_UP_2)

            # acquisition of the increase factor
            factor_accuracy_UP = input('\n\nEnter the INTEGER increase factor (>1): ')
            factor_accuracy_UP = int(factor_accuracy_UP)

            mesh_points_dummy = []
            p_new =[]

            # ***ATTENTION***
            # following this procedure 3 (...) different portions of the contour are considered:
            # (1) from 0 to the beginning of the more-resolved range where the contour is unchanged
            # (2) the more-resolved range in which the number of points is enhanced
            # (3) from the end of the more-resolved range to the end of the contour where it is unchanged
            # However, (1) and (3) are meaningless if, respectively: the more-resolved range starts from 0;
            # the more-resolved range ends in 0;

            # (1) from 0 to the beginning of the more-resolved range where the contour is unchanged
            for ii in range(0, where_accuracy_UP_1):

                p1 = mesh_points[ii]

                # storage of the first point of the ii-th segment
                mesh_points_dummy = mesh_points_dummy + [(p1)]

                #mesh_points_dummy = mesh_points_dummy + [(p2)]
                # NO NEED: it will be automatically inserted as the p1 at the step (ii + 1)-th


            # it is needed to separately considered the case in which where_accuracy_UP_2 is 0 or not:
            # indeed, in the former case, it is not possible to specify 
            # 'ii in range(where_accuracy_UP_1, where_accuracy_UP_2)' because where_accuracy_UP_2 = 0
            # and the range is not valid for a loop.
            # Hence, if where_accuracy_UP_2 = 0, the last step (from SIZE_MESH_POINTS - 1 to 0) is 
            # considered separately.
            end_point = 0
            flag = 0

            if where_accuracy_UP_2 != 0:
                end_point = where_accuracy_UP_2
            else:
                end_point = SIZE_MESH_POINTS - 1
                flag = 1

            
            # (2) the more-resolved range in which the number of points is enhanced
            for ii in range(where_accuracy_UP_1, end_point):

                p1 = mesh_points[ii]

                # storage of the first point of the ii-th segment
                mesh_points_dummy = mesh_points_dummy + [(p1)]

                if ii < (SIZE_MESH_POINTS - 1):
                    p2 = mesh_points[ii + 1]
                else :
                    p2 = mesh_points[0]  

                y_vector = p1.vector_to(p2)

                # new points creation: the vector between the extrema p1 and p2 is divided in factor_accuracy_UP parts
                # so that to identify the (factor_accuracy_UP - 1) inner points
                y_vector = y_vector / factor_accuracy_UP

                for jj in range(1, factor_accuracy_UP):

                    p_new = p1 + y_vector * jj
                    # storage of the jj-th new point of the ii-th segment
                    mesh_points_dummy = mesh_points_dummy + [(p_new)]

                # mesh_points_dummy = mesh_points_dummy + [(p2)]
                # NO NEED: it will be automatically inserted as the p1 at the step (ii + 1)-th

            # (3) from the end of the more-resolved range to the end of the contour where it is unchanged
            if flag == 0:
                for ii in range(where_accuracy_UP_2, SIZE_MESH_POINTS):

                    p1 = mesh_points[ii]

                    # storage of the first point of the ii-th segment
                    mesh_points_dummy = mesh_points_dummy + [(p1)]

                    #mesh_points_dummy = mesh_points_dummy + [(p2)]
                    # NO NEED: it will be automatically inserted as the p1 at the step (ii + 1)-th

            # (2) the more-resolved range in which the number of points is enhanced, LAST STEP, no need for (3)
            else:            
                p1 = mesh_points[-1]

                # storage of the first point of the ii-th segment
                mesh_points_dummy = mesh_points_dummy + [(p1)]

                p2 = mesh_points[0]

                y_vector = p1.vector_to(p2)
                y_vector = y_vector / factor_accuracy_UP

                for jj in range(1, factor_accuracy_UP):

                    p_new = p1 + y_vector * jj
                    # storage of the jj-th new point of the ii-th segment
                    mesh_points_dummy = mesh_points_dummy + [(p_new)]


            # the previous vector is overwritten
            mesh_points = mesh_points_dummy
            # to be updated! The new one is bigger than the old one
            SIZE_MESH_POINTS = int(np.size(mesh_points)) 
            
            mesh_points = from_Point2D_to_2xN_array(mesh_points)

            # plot of the new contour
            fig_1 = plt.figure(1)
            ax = fig_1.add_subplot(111)
            plt.plot(mesh_points[0,:], mesh_points[1,:], 'ko-')
            kk = 0
            for x,y in zip(mesh_points[0,:], mesh_points[1,:]):
                label = "{}".format(kk)
                kk = kk + 1
                plt.annotate(label, # this is the text
                             (x,y), # this is the point to label
                             textcoords="offset points", # how to position the text
                             xytext=(0,8), # distance from text to points (x,y)
                             ha='center') # horizontal alignment can be left, right or center
            plt.plot([mesh_points[0, -1], mesh_points[0, 0]], [mesh_points[1, -1], mesh_points[1, 0]], 'ko-')
            plt.title('New reference contour')
            plt.axis('equal')
            plt.show()

            mesh_points = from_2xN_array_to_Point2D(mesh_points)


        elif accuracy_UP != 'n' and accuracy_UP != 'N':
            print('\n\n***ERROR***\n\nInvalid input! The contour points are taken as reference, no increase in the accuracy!')


    #################################################################################################################
    #################################################################################################################

    #############
    # detectors #

    # the set of mesh_points is traveled CLOCKWISELY
    # it is chosen to have OUTWARD-POINTING normal unit vectors for the mesh points

    # distance between one mesh point and the corresponding one of the detectors' set
    ##########
    # resize #

    # the set of mesh_points is traveled CLOCKWISELY
    # it is chosen to have OUTWARD-POINTING normal unit vectors for the mesh points

    SIZE_CONTOUR_POINTS = SIZE_MESH_POINTS

    contour_points_new = mesh_points

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


    contour_points_resized = np.transpose(contour_points_resized)

    return contour_points_resized
