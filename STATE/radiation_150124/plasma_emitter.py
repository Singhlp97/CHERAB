
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
import os
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

# MMM 14.03.2022
import matplotlib.pyplot as plt
from import_EIRENE_data import load_eirene_data

# MMM 18.05.2022
from import_SOLPS_data import load_solps_data

def make_plasma_emitter(cfg = None, parent = None,
                        astra_simulation = None, solps_simulation = None):
    """
    Non-spectral emitter with the emissivity defined as SOLPSFunction2D.
    :param SOLPSMesh solps_mesh: SOLPS simulation mesh.
    :param SOLPSFunction2D radiation_function: Emissivity in W m-3.
    :param Node parent: parent node in the scenegraph, e.g. a World object.
    :param float step: Volume integration step in meters.
    :rtype: Primitive

    TODO mmm:
    - Check whether more efficient to have N 3D emitters or to have N 2D emitters
      then combined in 1 overall 3D emitter
    - Add possibility to choose among different type of radiation emission when
      Axisymmetric source is used
    """

    ####################################################
    ################## PRELIMINARIES ###################
    ####################################################

    RAD2DEG = 180 / np.pi

    use_ASTRA_emission    = cfg['plasma']['ASTRA']['use_ASTRA_emission']
    use_B2_emission       = cfg['plasma']['SOLPS']['use_B2_emission']
    use_EIRENE_emission   = cfg['plasma']['SOLPS']['use_EIRENE_emission']
    use_extra_emission    = cfg['plasma']['SOLPS']['use_extra_emission']

    # different weights [temporary]

    use_weighted_emission = cfg['plasma']['SOLPS']['weights']['use_weighted_emission']
    if use_weighted_emission is True:
        weight_main_emission  = cfg['plasma']['SOLPS']['weights']['weight_main_emission']
        weight_extra_emission = cfg['plasma']['SOLPS']['weights']['weight_extra_emission']

    #################################################
    ################## ABSORPTION ###################
    #################################################

    use_absorption_function = cfg['plasma']['absorption']['use_absorption_function']

    if use_absorption_function is True:
        (absorption_function_2d_astra,
         absorption_function_2d_b2,
         absorption_function_2d_eirene,
         absorption_function_2d_extra) = CreateAbsorptionFunction(cfg = cfg,
                                         solps_simulation = solps_simulation)
        absorption_function_3d = 0.0
    else:
        absorption_function_2d = None
        absorption_function_3d = None

    #####################################################
    ################## EXTRA B2 EMISSION ################
    #####################################################

    if use_extra_emission is True:
        solps_extra = load_solps_data(cfg = cfg, extra = True)

    #####################################################
    ################## INSIDE B2 MESH ###################
    #####################################################

    inside_mesh = solps_simulation._inside_mesh

    # add inside_mesh from solps_extra where original inside_mesh is zero
    # to avoid overlapping (i.e. resulting inside_mesh = 2)

    if use_extra_emission is True:
        inside_mesh += solps_extra._inside_mesh * (1.0 - inside_mesh)

    #####################################################
    ################## OUTSIDE B2 MESH ##################
    #####################################################

    # EIRENE data employed are those of original solps_simulation
    # NOT those from solps_extra

    if use_EIRENE_emission is True:

        eirene_simulation = load_eirene_data(cfg = cfg)
        outside_B2_mesh  = 1.0 - solps_simulation._inside_mesh
        if use_extra_emission is True: outside_B2_mesh *= (1.0 - solps_extra._inside_mesh)

    #####################################################
    ################## TYPE OF MAPPER ###################
    #####################################################

    # CAUTION.
    #
    # - ASTRA & EIRENE emission from neutrals outside B2 mesh exist all along the toroidal direction
    # - ASTRA & EIRENE emission will be sampled with step_max because:
    #   - Un-structured EIRENE mesh => knowledge about neighboring cells but COMPLICATED
    #   - Mild gradients expected anyway

    ##################################

    if use_B2_emission is True:

        if cfg['plasma']['type_radiation'] == "halpha_total_radiation":
            solps_radiation_f2d = solps_simulation.halpha_total_radiation_f2d

        elif cfg['plasma']['type_radiation'] == "halpha_mol_radiation":
            solps_radiation_f2d = solps_simulation.halpha_mol_radiation_f2d

        elif cfg['plasma']['type_radiation'] == "total_radiation":
            solps_radiation_f2d = solps_simulation.total_radiation_f2d

    ##################################

    try:    astra_radiation_f2d  = astra_simulation.Ha_emission_f2d
    except: astra_radiation_f2d  = 0.0 * solps_radiation_f2d
    try:    eirene_radiation_f2d = eirene_simulation.Ha_emission_f2d * outside_B2_mesh
    except: eirene_radiation_f2d = 0.0 * solps_radiation_f2d
    try:    solps_extra_f2d      = solps_extra.halpha_total_radiation_f2d * (1.0 - solps_simulation._inside_mesh)
    except: solps_extra_f2d      = 0.0 * solps_radiation_f2d

    if cfg['plasma']['Mapper'] == "DiscreteToroidalMapper":
        
        periodicity = 2 * np.pi / int(cfg['limiters']['number'])
        where_non_zero = cfg['limiters']['angular_width_deg'] / RAD2DEG

        try:
            # MMM 2022-07-27: ASTRA emission made discrete too
            radiation_function_3d = use_B2_emission * \
                                     DiscreteToroidalMapper(solps_radiation_f2d,
                                                            periodicity = periodicity,
                                                            where_non_zero = where_non_zero) + \
                                     use_ASTRA_emission * \
                                     DiscreteToroidalMapper(astra_radiation_f2d,
                                                            periodicity = periodicity,
                                                            where_non_zero = where_non_zero)
        except: pass

        # radiation_function_3d += AxisymmetricMapper(use_ASTRA_emission  * astra_radiation_f2d  + \
        #                                             use_EIRENE_emission * eirene_radiation_f2d + \
        #                                             use_extra_emission  * solps_extra_f2d)

        radiation_function_3d += AxisymmetricMapper(use_EIRENE_emission * eirene_radiation_f2d + \
                                                    use_extra_emission  * solps_extra_f2d)

        if use_absorption_function is True:

            # CAUTION.
            #
            # hp: step_max to sample ASTRA, otherwise absorption_function_2d must be revisited

            try: absorption_function_3d = use_B2_emission * \
                                          DiscreteToroidalMapper(absorption_function_2d_b2,
                                                                 periodicity = periodicity,
                                                                 where_non_zero = where_non_zero) 
            except: pass

            try:                 
                absorption_function_3d += AxisymmetricMapper(use_ASTRA_emission  * absorption_function_2d_astra + \
                                                             use_EIRENE_emission * absorption_function_2d_eirene + 
                                                             use_extra_emission  * absorption_function_2d_extra)
            except: pass

    elif cfg['plasma']['Mapper'] == "AxisymmetricMapper":
        
        radiation_function_3d = AxisymmetricMapper(use_ASTRA_emission  * astra_radiation_f2d  +
                                                   use_B2_emission     * solps_radiation_f2d  +
                                                   use_EIRENE_emission * eirene_radiation_f2d +
                                                   use_extra_emission  * solps_extra_f2d)

        if use_absorption_function is True:

            absorption_function_3d = AxisymmetricMapper(use_ASTRA_emission  * absorption_function_2d_astra +
                                                        use_B2_emission     * absorption_function_2d_b2    +
                                                        use_EIRENE_emission * absorption_function_2d_eirene + 
                                                        use_extra_emission  * absorption_function_2d_extra)

    #######################################################
    ################## TYPE OF SAMPLING ###################
    #######################################################

    use_step_function = cfg['raytracing']['sampling']['use_step_function']

    if use_step_function is True:
        step_function_3d = CreateStepFunction(cfg = cfg,
                                              solps_simulation = solps_simulation)
        if use_extra_emission is True:
            step_extra_3d = CreateStepFunction(cfg = cfg,
                                               solps_simulation = solps_simulation,
                                               extra = True)
            step_function_3d += step_extra_3d
    else:
        step_function_3d = None

    #########################################################
    ################## RADIATIVE EMISSION ###################
    #########################################################

    emitter = RadiationFunction(radiation_function =      radiation_function_3d,

                                use_step_function =       use_step_function,
                                step_function_3d =        step_function_3d,
                                step_max =                cfg['raytracing']['sampling']['step_max'],

                                use_scattering_function = False, # test
                                scattering_function_3d =  None,  # test: scattering
                                collisions_max =          1,     # test: 100 by default

                                use_absorption_function = use_absorption_function,
                                absorption_function_3d =  absorption_function_3d,

                                step = cfg['raytracing']['sampling']['step_uniform'])

    material = VolumeTransform(emitter, transform = translate(0, 0, 0))

    ##########################################################
    ################## ENCLOSING PRIMITIVE ###################
    ##########################################################

    solps_mesh = solps_simulation.mesh

    max_r = solps_mesh.mesh_extent['maxr']
    min_r = solps_mesh.mesh_extent['minr']
    max_z = solps_mesh.mesh_extent['maxz']
    min_z = solps_mesh.mesh_extent['minz']

    if use_extra_emission is True:
        extra_mesh = solps_extra.mesh
        max_r = np.max([max_r, extra_mesh.mesh_extent['maxr']])
        min_r = np.min([min_r, extra_mesh.mesh_extent['minr']])
        max_z = np.max([max_z, extra_mesh.mesh_extent['maxz']])
        min_z = np.min([min_z, extra_mesh.mesh_extent['minz']])

    if use_EIRENE_emission is True:
        eirene_mesh = eirene_simulation.mesh
        max_r = np.max([max_r, eirene_mesh['mesh_extent']['maxr']])
        min_r = np.min([min_r, eirene_mesh['mesh_extent']['minr']])
        max_z = np.max([max_z, eirene_mesh['mesh_extent']['maxz']])
        min_z = np.min([min_z, eirene_mesh['mesh_extent']['minz']])

    # padding

    padding = 1E-04
    max_r += padding
    min_r -= padding
    max_z += padding
    min_z -= padding

    # hollow cylinder

    height = max_z - min_z
    hollow_cylinder = Subtract(Cylinder(max_r, height),
                               Cylinder(min_r, height),
                               transform = translate(0, 0, +min_z))

    ######################################################
    ################## CREATING SOURCE ###################
    ######################################################   

    if cfg['plasma']['Mapper'] == "DiscreteToroidalMapper":
        align = - 0.5 * RAD2DEG * (periodicity - where_non_zero)
    else:
        align = 0.0

    plasma_volume = hollow_cylinder
    plasma_volume.parent = parent
    material.transform  = translate(0, 0, -min_z) * rotate_z(align)
    plasma_volume.material = material

    ######################################################
    ################## PLOTTING SOURCE ###################
    ######################################################

    if cfg['plotting']['plot_total_emission'] is True:

        num = 400
        ri = np.linspace(min_r, max_r, num)
        zi = np.linspace(min_z, max_z, num)

        # mind imshow convention

        emission = np.zeros((np.size(zi), np.size(ri)))
        emission_f3d = radiation_function_3d

        for i in range(emission.shape[0]):
            for j in range(emission.shape[1]):
                emission[i,j] = emission_f3d(ri[j], 0.0, zi[i])

        fig, ax = plt.subplots()
        c = ax.pcolormesh(ri, zi, np.log10(emission + emission.min() * 0.2), cmap = 'jet')
        fig.colorbar(c, ax = ax)
        ax.set_title('Radiative emission')
        # ax.set_xlim(left = 0.1, right = 1.1)
        ax.axis('equal')

        if cfg['plotting']['save_figures'] is True:
            if os.path.exists(os.path.join(cfg['output_directory_extended'], 'emission.png')) is True:
                plt.savefig(os.path.join(cfg['output_directory_extended'], 'emission_extra.png'))
            else:
                plt.savefig(os.path.join(cfg['output_directory_extended'], 'emission.png'))

        #######################

    return emitter
