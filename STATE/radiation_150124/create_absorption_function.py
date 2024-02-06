
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

# MMM
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io

from import_ASTRA_data  import load_astra_data
from import_SOLPS_data  import load_solps_data
from import_EIRENE_data import load_eirene_data

from quantum_data import FractionExcitedHydrogenPlasma, FractionExcitedHydrogenGas

from homemade_distributions import make_homemade_emission

from cherab.solps.solps_2d_functions import SOLPSFunction2D
from cherab.core.math import Interpolate2DLinear
from cherab.tools.equilibrium import import_eqdsk

from accessories import plot_2d, make_chamber_wall

"""
TODO mmm:
- interpolation to get absorption_function_2d_core to be done in ASTRASimulation
  object, likewise in EIRENESimulation
- clean this mess...
"""

def CreateAbsorptionFunction(cfg = None, solps_simulation = None):

    ################

    use_ASTRA_emission  = cfg['plasma']['ASTRA']['use_ASTRA_emission']
    use_B2_emission  = cfg['plasma']['SOLPS']['use_B2_emission']
    use_EIRENE_emission = cfg['plasma']['SOLPS']['use_EIRENE_emission']
    use_extra_emission = cfg['plasma']['SOLPS']['use_extra_emission']
    use_homemade_emission = cfg['plasma']['homemade']['use_homemade_emission']

    # MMM: different weights
    use_weighted_emission = cfg['plasma']['SOLPS']['weights']['use_weighted_emission']
    if use_weighted_emission is True:
        # only weight of main emission by assumption (ASTRA from main emission)
        weight_main  = cfg['plasma']['SOLPS']['weights']['weight_main_emission']
        weight_extra = cfg['plasma']['SOLPS']['weights']['weight_extra_emission']

    plot_absorption_function = cfg['plotting']['plot_absorption_function']

    absorption_function_2d_core     = 0.0
    absorption_function_2d_sol      = 0.0
    absorption_function_2d_eir      = 0.0
    absorption_function_2d_extra    = 0.0

    ###########################################################################################################################################
    sigma_Ha =  cfg['plasma']['absorption']['absorption_cross_section_H23'] # 1.195E-14 [m^2]: from absorption/Ha_absorption_cross-section.py #
    ###########################################################################################################################################

    if use_extra_emission is True:
        solps_extra = load_solps_data(cfg = cfg, extra = True)

    ################

    if use_ASTRA_emission is True:

        astra_simulation = load_astra_data(cfg = cfg)

        ng_core = astra_simulation.neutral_density_raw
        ne_core = astra_simulation.electron_density_raw
        Te_core = astra_simulation.electron_temperature_raw
        fexct_core = FractionExcitedHydrogenPlasma(n_m3 = ne_core, T_eV = Te_core)
        fexct_core[np.isnan(fexct_core)] = 0
        ng_2_core = ng_core * fexct_core
        absorption_function_2d_raw_core = sigma_Ha * ng_2_core
        f_excited_H_core_2d_raw = fexct_core

        if use_weighted_emission is True: absorption_function_2d_raw_core *= weight_main

    ###################

    if use_EIRENE_emission is True:

        eirene_simulation = load_eirene_data(cfg = cfg)

        # H2 and H2+ currently neglected - but possibly NOT to be includede at all

        ng_eir = eirene_simulation.pdena
        Tg_eir = eirene_simulation.tdena
        fexct_eir = FractionExcitedHydrogenGas(n_quantum = 2, T_eV = Tg_eir) # H*(2) fractional abundance
        ng_2_eir = ng_eir * fexct_eir
        absorption_function_2d_raw_eir = sigma_Ha * ng_2_eir

        eirene_simulation.SigmaAbs = sigma_Ha * ng_2_eir
        absorption_function_2d_eir = eirene_simulation.SigmaAbs_f2d * (1.0 - solps_simulation._inside_mesh) # * 1E+04

        if use_extra_emission is True:
            absorption_function_2d_eir *= (1.0 - solps_extra._inside_mesh)

        eirene_simulation.f_exct = fexct_eir
        f_excited_H_eir_2d_raw = fexct_eir
        f_excited_H_eir_2d = eirene_simulation.f_exct_f2d * (1.0 - solps_simulation._inside_mesh)

        if use_weighted_emission is True: absorption_function_2d_eir *= weight_main

    ###################

    if use_B2_emission is True:

        # sol
        # raw data on B2.5 (nx, ny) grid

        fexct = FractionExcitedHydrogenPlasma(n_m3 = solps_simulation.electron_density,
                                              T_eV = solps_simulation.electron_temperature)
        ng_sol = solps_simulation.species_density[0,:,:] # mind CHERAB index convention
        ng_2_sol = ng_sol * fexct
        absorption_function_2d_raw_sol = sigma_Ha * ng_2_sol
        f_excited_H_sol_2d_raw = fexct

        if use_weighted_emission is True: absorption_function_2d_raw_sol *= weight_main

        if use_extra_emission is True:

            solps_extra = load_solps_data(cfg = cfg, extra = True)
            fexct = FractionExcitedHydrogenPlasma(n_m3 = solps_extra.electron_density,
                                                  T_eV = solps_extra.electron_temperature)
            ng_extra = solps_extra.species_density[0,:,:] # mind CHERAB index convention
            ng_2_extra = ng_extra * fexct
            absorption_function_2d_raw_extra = sigma_Ha * ng_2_extra
            f_excited_H_extra_2d_raw = fexct

            if use_weighted_emission is True: absorption_function_2d_raw_extra *= weight_extra

    ###################

    if use_homemade_emission is True:

        ((ng_2_f2d, ng_2_raw_2d), _, _) = make_homemade_emission(cfg = cfg, verbose = False)
        absorption_function_2d = ng_2_f2d * sigma_Ha

    #############################
    # sn_max for Delta-tracking #
    #############################

    # Maximum can be computed from raw data:
    # linear interpolation does NOT lead to any value above raw maximum

    if use_homemade_emission is True:
        sn_max = ng_2_raw_2d.max() * sigma_Ha
    else:
        raise ValueError('To be implemented!')
        # .max() method does not work if any is 0.0 (float)
        sn_max = np.max([absorption_function_2d_raw_core.max(),
                         absorption_function_2d_raw_sol.max(),
                         absorption_function_2d_raw_eir.max()])

    ####################
    # 2D interpolation #
    ####################

    if use_homemade_emission is False:
        
        # core
        if use_ASTRA_emission is True:
            r = astra_simulation.equilibrium.r_data
            z = astra_simulation.equilibrium.z_data
            a = absorption_function_2d_raw_core.T
            f = f_excited_H_core_2d_raw.T
            outside_SOLPS = (1.0 - solps_simulation._inside_mesh) * astra_simulation.equilibrium.inside_lcfs
            absorption_function_2d_core = Interpolate2DLinear(r, z, a, extrapolate = True) * outside_SOLPS
            f_excited_H_core_2d = Interpolate2DLinear(r, z, f, extrapolate = True) * outside_SOLPS

        # sol
        absorption_function_2d_sol = SOLPSFunction2D.instance(solps_simulation._inside_mesh,
                                                              absorption_function_2d_raw_sol)
        f_excited_H_sol_2d = SOLPSFunction2D.instance(solps_simulation._inside_mesh,
                                                              f_excited_H_sol_2d_raw)

        if use_extra_emission is True:
            absorption_function_2d_extra = SOLPSFunction2D.instance(solps_extra._inside_mesh,
                                                                    absorption_function_2d_raw_extra)
            f_excited_H_extra_2d = SOLPSFunction2D.instance(solps_extra._inside_mesh,
                                                            f_excited_H_extra_2d_raw)        
            absorption_function_2d_extra *= (1.0 - solps_simulation._inside_mesh)
            f_excited_H_extra_2d         *= (1.0 - solps_simulation._inside_mesh)

        absorption_function_2d = absorption_function_2d_sol
        f_excited_2d = f_excited_H_sol_2d

        if use_ASTRA_emission is True:
            absorption_function_2d += absorption_function_2d_core
            f_excited_2d += f_excited_H_core_2d
        if use_EIRENE_emission is True:
            absorption_function_2d += absorption_function_2d_eir
            f_excited_2d += f_excited_H_eir_2d
        if use_extra_emission is True:
            absorption_function_2d += absorption_function_2d_extra
            f_excited_2d += f_excited_H_extra_2d

    ########################

    if plot_absorption_function is True:

        if use_homemade_emission is False:
            mesh = solps_simulation.mesh
            minr = 0.2 # solps_simulation.mesh.mesh_extent['minr']
            maxr = 1.1 # solps_simulation.mesh.mesh_extent['maxr']
            minz = -0.9 #solps_simulation.mesh.mesh_extent['minz']
            maxz = 0.9 #solps_simulation.mesh.mesh_extent['maxz']
        else:
            location = os.path.join(cfg['baserun'],
                                    'input',
                                    cfg['run'],
                                    cfg['plasma']['ASTRA']['ASTRA_directory'])
            eqFile = cfg['plasma']['ASTRA']['geqdsk_name']
            eqTime = cfg['plasma']['ASTRA']['geqdsk_time']
            equilibrium = import_eqdsk(os.path.join(location, eqFile + '_' + eqTime + '.geqdsk'))

        plot_2d(ri = equilibrium.r_data, zi = equilibrium.z_data, f2d = absorption_function_2d, scale = "log", vmin = 1E-04, eq = equilibrium, title = r'Absorption macroscopic cross-section $\Sigma=\sigma\cdot n_{2}^{*}$', cmap_label = '$log_{10} \\; [m^{-1}]$')

        if cfg['plotting']['save_figures'] is True:
            plt.savefig(os.path.join(cfg['output_directory_extended'], 'abs_Sigma.png'))

        if use_homemade_emission is False:
            plot_2d(ri = equilibrium.r_data, zi = equilibrium.z_data, f2d = f_excited, scale = "log", eq = None, title = r'Fractional abundance $H(2)/H_{tot}$')
            if cfg['plotting']['save_figures'] is True:
                plt.savefig(os.path.join(cfg['output_directory_extended'], 'abs_fexct.png'))

        #######################

        #plt.show()

    if use_homemade_emission is False:
        return (absorption_function_2d_core, absorption_function_2d_sol,
                absorption_function_2d_eir, absorption_function_2d_extra)
    else:
        return (absorption_function_2d, sn_max)

############################################################################################
############################################################################################
############################################################################################