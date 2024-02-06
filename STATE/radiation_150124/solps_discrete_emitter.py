
# Copyright 2016-2021 Euratom
# Copyright 2016-2021 United Kingdom Atomic Energy Authority
# Copyright 2016-2021 Centro de Investigaciones Energéticas, Medioambientales y Tecnológicas
#
# Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the
# European Commission - subsequent versions of the EUPL (the "Licence");
# You may not use this work except in compliance with the Licence.
# You may obtain a copy of the Licence at:
#
# https://joinup.ec.europa.eu/software/page/eupl5
#
# Unless required by applicable law or agreed to in writing, software distributed
# under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied.
#
# See the Licence for the specific language governing permissions and limitations
# under the Licence.

from raysect.core import translate
from raysect.optical.material import VolumeTransform
from raysect.primitive import Cylinder, Subtract

from cherab.core.math.mappers import AxisymmetricMapper, DiscreteToroidalMapper
from cherab.tools.emitters import RadiationFunction

from cherab.solps import SOLPSMesh, SOLPSFunction2D

# MMM
import numpy as np                
from raysect.core import Point3D, rotate_z
from raysect.primitive import Box

# MMM 24.01.2022
from raysect.core import Vector3D
from raysect.primitive import Intersect

# MMM 27.01.2022
from create_step_function import CreateStepFunction

# MMM 29.01.2022
from create_absorption_function import CreateAbsorptionFunction

def make_solps_discrete_emitter(solps_mesh = None,
                                solps_simulation = None, radiation_function = None,
                                parent = None, step = 0.01, cfg = None):
    """
    Non-spectral emitter with the emissivity defined as SOLPSFunction2D.
    :param SOLPSMesh solps_mesh: SOLPS simulation mesh.
    :param SOLPSFunction2D radiation_function: Emissivity in W m-3.
    :param Node parent: parent node in the scenegraph, e.g. a World object.
    :param float step: Volume integration step in meters.
    :rtype: Primitive
    """

    ####################################################
    ################## PRELIMINARIES ###################
    ####################################################

    if not isinstance(solps_mesh, SOLPSMesh):
        raise TypeError('Argument solps_mesh must be a SOLPSMesh instance.')
    if not isinstance(radiation_function, SOLPSFunction2D):
        raise TypeError('Argument radiation_function must be a SOLPSFunction2D instance.')

    RAD2DEG = 180 / np.pi

    #################################################
    ################## ABSORPTION ###################
    #################################################

    use_absorption_function = cfg['plasma']['absorption']['use_absorption_function']

    print()
    print("######################")
    print("use_absorption_function = {}".format(use_absorption_function))

    if use_absorption_function == 1:
        absorption_function_2d = CreateAbsorptionFunction(cfg = cfg,
                                                          solps_simulation = solps_simulation)
    else:
        absorption_function_2d = None
        absorption_function_3d = None

    #####################################################
    ################## TYPE OF MAPPER ###################
    #####################################################

    if cfg['plasma']['Mapper'] == "DiscreteToroidalMapper":
        
        periodicity = 2 * np.pi / int(cfg['limiters']['number'])
        where_non_zero = cfg['limiters']['angular_width_deg'] / RAD2DEG

        radiation_function_3d = DiscreteToroidalMapper(radiation_function,
                                                       periodicity = periodicity,
                                                       where_non_zero = where_non_zero)            

        if use_absorption_function == 1:
            absorption_function_3d = DiscreteToroidalMapper(absorption_function_2d,
                                                            periodicity = periodicity,
                                                            where_non_zero = where_non_zero)            

    elif cfg['plasma']['Mapper'] == "AxisymmetricMapper":

        radiation_function_3d = AxisymmetricMapper(radiation_function)

        if use_absorption_function == 1:
            absorption_function_3d = AxisymmetricMapper(absorption_function_2d)

    #######################################################
    ################## TYPE OF SAMPLING ###################
    #######################################################

    use_step_function = cfg['plasma']['sampling']['use_step_function']

    print()
    print("######################")
    print("use_step_function = {}".format(use_step_function))

    if use_step_function == 1:
        step_function_3d = CreateStepFunction(cfg = cfg,
                                              solps_simulation = solps_simulation,
                                              emission_3d = radiation_function_3d)
    else:
        step_function_3d = None

    #########################################################
    ################## RADIATIVE EMISSION ###################
    #########################################################

    emitter = RadiationFunction(radiation_function =      radiation_function_3d,

                                use_step_function =       use_step_function,
                                step_function_3d =        step_function_3d,
                                step_max =                cfg['plasma']['sampling']['step_max'],

                                use_scattering_function = 0,     # test
                                scattering_function_3d =  None,  # test: scattering ~ absorption
                                collisions_max =          1,     # test: 100 by default

                                use_absorption_function = use_absorption_function,
                                absorption_function_3d =  absorption_function_3d,

                                step = step)

    material = VolumeTransform(emitter, transform = translate(0, 0, 0))

    ##########################################################
    ################## ENCLOSING PRIMITIVE ###################
    ##########################################################

    outer_radius = solps_mesh.mesh_extent['maxr']
    inner_radius = solps_mesh.mesh_extent['minr']
    height = solps_mesh.mesh_extent['maxz'] - solps_mesh.mesh_extent['minz']
    lower_z = solps_mesh.mesh_extent['minz']
    min_z = solps_mesh.mesh_extent['minz']
    max_z = solps_mesh.mesh_extent['maxz']

    # hollow cylinder
    hollow_cylinder = Subtract(Cylinder(outer_radius, height),
                               Cylinder(inner_radius, height),
                               transform = translate(0, 0, +lower_z))

    ######################################################
    ################## CREATING SOURCE ###################
    ######################################################   

    print()
    print("######################")
    print("Mapper = {}".format(cfg['plasma']['Mapper']))
    print()

    if cfg['plasma']['Mapper'] == "DiscreteToroidalMapper":

        align = - 0.5 * RAD2DEG * (periodicity - where_non_zero)

        plasma_volume = hollow_cylinder
        plasma_volume.parent = parent
        material.transform  = translate(0, 0, -lower_z) * rotate_z(align)
        plasma_volume.material = material

    elif cfg['plasma']['Mapper'] == "AxisymmetricMapper":
        

        # conventional toroidal plasma emitter # MMM
        # plasma_volume = Subtract(Cylinder(outer_radius, height), Cylinder(inner_radius, height),
        #                         material=material, parent=parent, transform=translate(0, 0, lower_z))

        ###########################################################
        # MMM whatever comes next
        #
        # CHERAB sees an emitter if World contains at least one
        # (material + enclosing primitive), i.e. whatever outside
        # a primitive boundary is NOT seen as a source.
        #
        # Leading idea:
        #
        # - create a toroidally symmetric source of radiation
        # - create discrete slices of the usual hollow cylinder
        #   (toroidally symmetric object)
        # 
        # => toroidally symmetric source 
        #    +
        #    discrete primitives
        #    =
        #    discrete source 
        ###########################################################

        # BOX
        #
        # Create box to slice hollow cylinder:
        # - one side through the origin 
        # - linear size (L) > max(hollow_cylinder height, hollow_cylinder diameter)
        # - translated so as mid-height at R-phi plane (Z=0)
        # - lower_point: initially in origin , then translated to actual location
        # - upper_point: initially in (1,1,1), then translated to actual location
        # - rigid body translation <=> same translation for both lower_ and upper_point
        # - emission exist only on limiters (i.e. where plasma touches walls)
        # - angular_width_sector = angular_width_lim + angular_width_nonlim
        #   and
        #   angular_width_sector_deg = 360 / num_limiters

        deg_to_rad = np.pi / 180                                     # conversion factor

        num_limiters = cfg['limiters']['number']                     # [ - ] 8 limiters = 8 non-limiters
        angular_width_lim_deg = cfg['limiters']['angular_width_deg'] # [deg] width of the limiters ????

        angular_width_lim_rad = angular_width_lim_deg * deg_to_rad   # [rad] 
        angular_width_sector_deg = 360 / num_limiters                # [deg]

        L = np.max([height, 2 * outer_radius]) * 1E+01 # [m]: > max(hollow_cylinder height, hollow_cylinder diameter)

        # box limiting points
        to_right_location = - L * Vector3D(0.5, 1.0, 0.5)
        lower_point = Point3D(0, 0, 0) + to_right_location
        upper_point = Point3D(L, L, L) + to_right_location

        # 1st box: 180 degrees rotation because I like it 
        box_1 = Box(lower_point, upper_point, transform = rotate_z(180))
        
        # subtract box from hollow_cylinder (sliced in half) => half_hollow_cylinder
        half_hollow_cylinder = Subtract(hollow_cylinder, box_1)
        
        # 2nd box: must be properly aligned ~ angular_width_lim_deg
        box_2 = Box(lower_point, upper_point, transform = rotate_z(180 - angular_width_lim_deg))
        
        # for loop to create num_limiters plasma volumes
        for i in range(num_limiters):
            # Intersection (not subtraction...) between half_hollow_cylinder and box_2 gives the slice
            plasma_volume = Intersect(half_hollow_cylinder, box_2, 
                                      parent = parent, material = material,
                                      transform = rotate_z(i * angular_width_sector_deg))

    # no need to return anything: plasma_volumes already in world

    return emitter #plasma_volume #step_function_3d
