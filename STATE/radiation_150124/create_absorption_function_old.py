
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

from cherab.solps.solps_2d_functions import SOLPSFunction2D
from cherab.core.math import Interpolate2DLinear

"""
TODO mmm:
- interpolation to get absorption_function_2d_core to be done in ASTRASimulation
  object, likewise in EIRENESimulation
- clean this mess...
"""

def CreateAbsorptionFunction(cfg = None, solps_simulation = None):

    ################

    use_ASTRA_emission  = cfg['plasma']['ASTRA']['use_ASTRA_emission']
    use_EIRENE_emission = cfg['plasma']['SOLPS']['use_EIRENE_emission']

    plot_absorption_function = cfg['plotting']['plot_absorption_function']

    ######################
    sigma_Ha = 1.195E-14 # [m^2]: from absorption/Ha_absorption_cross-section.py
    ######################

    ################

    if use_ASTRA_emission is True:

        astra_simulation = load_astra_data(cfg = cfg)

        nH0_core = astra_simulation.neutral_density_raw
        ne_core = astra_simulation.electron_density_raw
        Te_core = astra_simulation.electron_temperature_raw
        fexct_core = FractionExcitedHydrogenPlasma(n_m3 = ne_core, T_eV = Te_core)
        fexct_core[np.isnan(fexct_core)] = 0
        nH0_2_core = nH0_core * fexct_core
        absorption_function_2d_raw_core = sigma_Ha * nH0_2_core
        f_excited_H_core_2d_raw = fexct_core

    if use_EIRENE_emission is True:

        eirene_simulation = load_eirene_data(cfg = cfg)

        # H2 and H2+ currently neglected

        nH0_eir = eirene_simulation.pdena
        TH0_eir = eirene_simulation.tdena
        fexct_eir = FractionExcitedHydrogenGas(n_quantum = 2, T_eV = TH0_eir) # H*(2) fractional abundance
        nH0_2_eir = nH0_eir * fexct_eir
        absorption_function_2d_raw_eir = sigma_Ha * nH0_2_eir

        eirene_simulation.SigmaAbs = sigma_Ha * nH0_2_eir
        absorption_function_2d_eir = eirene_simulation.SigmaAbs_f2d * (1.0 - solps_simulation._inside_mesh) # * 1E+04

        eirene_simulation.f_exct = fexct_eir
        f_excited_H_eir_2d_raw = fexct_eir
        f_excited_H_eir_2d = eirene_simulation.f_exct_f2d * (1.0 - solps_simulation._inside_mesh)

    # sol
    # raw data on B2.5 (nx, ny) grid

    fexct = FractionExcitedHydrogenPlasma(n_m3 = solps_simulation.electron_density,
                                          T_eV = solps_simulation.electron_temperature)
    nH0_sol = solps_simulation.species_density[0,:,:] # mind CHERAB index convention
    nH0_2_sol = nH0_sol * fexct
    absorption_function_2d_raw_sol = sigma_Ha * nH0_2_sol
    f_excited_H_sol_2d_raw = fexct

    ####################
    # 2D interpolation #
    ####################
    
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

    absorption_function_2d = absorption_function_2d_sol
    f_excited_2d = f_excited_H_sol_2d

    if use_ASTRA_emission is True:
        absorption_function_2d += absorption_function_2d_core
        f_excited_2d += f_excited_H_core_2d
    if use_EIRENE_emission is True:
        absorption_function_2d += absorption_function_2d_eir
        f_excited_2d += f_excited_H_eir_2d

    ########################

    if plot_absorption_function is True:

        mesh = solps_simulation.mesh

        num = 500
        minr = 0.2 # solps_simulation.mesh.mesh_extent['minr']
        maxr = 1.1 # solps_simulation.mesh.mesh_extent['maxr']
        minz = -0.9 #solps_simulation.mesh.mesh_extent['minz']
        maxz = 0.9 #solps_simulation.mesh.mesh_extent['maxz']
        ri = np.linspace(minr, maxr, num)
        zi = np.linspace(minz, maxz, num)

        if cfg['plasma']['type_radiation'] == "halpha_total_radiation":
            solps_radiation_f2d = solps_simulation.halpha_total_radiation_f2d
        elif cfg['plasma']['type_radiation'] == "halpha_mol_radiation":
            solps_radiation_f2d = solps_simulation.halpha_mol_radiation_f2d
        elif cfg['plasma']['type_radiation'] == "total_radiation":
            solps_radiation_f2d = solps_simulation.total_radiation_f2d

        # mind imshow convention
        absorption = np.zeros((np.size(zi), np.size(ri)))
        f_excited = np.zeros((np.size(zi), np.size(ri)))
        emission = np.zeros((np.size(zi), np.size(ri)))
        emission_2d = solps_radiation_f2d
        if use_ASTRA_emission is True:
            emission_2d += astra_simulation.Ha_emission_f2d * outside_SOLPS
        if use_EIRENE_emission is True:
            emission_2d += eirene_simulation.Ha_emission_f2d * (1.0 - solps_simulation._inside_mesh)
        for i in range(absorption.shape[0]):
            for j in range(absorption.shape[1]):
                absorption[i,j] = absorption_function_2d(ri[j], zi[i])
                f_excited[i,j] = f_excited_2d(ri[j], zi[i])
                emission[i,j] = emission_2d(ri[j], zi[i])
        absorption[absorption < absorption.max() / 1E+05] = 0
        f_excited[f_excited < f_excited.max() / 1E+05] = 0
    
        fig, ax = plt.subplots()
        c = ax.pcolormesh(ri, zi, np.log10(absorption + absorption.min() * 0.2))
        fig.colorbar(c, ax = ax)
        ax.set_title(r'Absorption macroscopic cross-section $\Sigma=\sigma\cdot n_{2}^{*}$')
        ax.axis('equal')

        if cfg['plotting']['save_figures'] is True:
            plt.savefig(os.path.join(cfg['output_directory_extended'], 'abs_Sigma.png'))

        fig, ax = plt.subplots()
        c = ax.pcolormesh(ri, zi, np.log10(f_excited + f_excited.min() * 0.2))
        fig.colorbar(c, ax = ax, cmap = 'jet')
        ax.set_title(r'Fractional abundance $H(2)/H_{tot}$')
        ax.axis('equal')

        if cfg['plotting']['save_figures'] is True:
            plt.savefig(os.path.join(cfg['output_directory_extended'], 'abs_fexct.png'))

        #######################

        plt.show()
    
    return absorption_function_2d

############################################################################################
############################################################################################
############################################################################################