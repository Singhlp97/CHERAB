
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
import pickle      
import matplotlib.pyplot as plt  
from matplotlib.collections import PolyCollection    

from import_SOLPS_data import load_solps_data, load_solps_mesh_padded
from cherab.solps import SOLPSSimulation
from cherab.solps.solps_2d_functions import SOLPSFunction2D

# MMM 30.01.2022
from differential_operations import ComputeGradientNorm, ComputeHessianNorm2
from cherab.core.math import Interpolate1DLinear, Interpolate2DLinear
from cherab.core.math.mappers import AxisymmetricMapper

# MMM 18.05.2022
from import_SOLPS_data import load_solps_data
from solps_python_scripts.read_b2fgmtry import read_b2fgmtry

from convert import homemade_to_solps_format
from accessories import plot_2d, make_chamber_wall

####################################################################
####################################################################

def CreateStepFunction(cfg = None, solps_simulation = None, extra = False):

    step_max = cfg['raytracing']['sampling']['step_max']
    step_min = cfg['raytracing']['sampling']['step_min']
    alpha    = cfg['raytracing']['sampling']['alpha']
    beta     = cfg['raytracing']['sampling']['beta']

    plot_step_function = cfg['plotting']['plot_step_function']

    use_homemade_emission = cfg['plasma']['homemade']['use_homemade_emission']

    if use_homemade_emission is False:
        SOLPS_directory = os.path.join(cfg['baserun'], 'input', cfg['run'], cfg['plasma']['SOLPS']['SOLPS_directory'])
        b2fgmtry = read_b2fgmtry(where = SOLPS_directory, verbose = False, save = True)
    elif use_homemade_emission is True:
        (b2fgmtry, solps_simulation) = homemade_to_solps_format(cfg = cfg)

    ################

    # CAUTION.
    #
    # ASTRA - core plasma: step = step_max = constant is assumed
    
    # solps_simulation = load_solps_data(cfg = cfg)
    
    # compute actual maximum from SOLPS simulation raw data

    if cfg['plasma']['type_radiation'] == "halpha_total_radiation":
        if use_homemade_emission is False:
            emission_2d_raw = solps_simulation.halpha_total_radiation
            if extra is True:
                solps_extra = load_solps_data(cfg = cfg, extra = True)
                # CAUTION.
                #
                # - ideally: inside_mesh = (1.0 - solps_simulation._inside_mesh) * solps_extra._inside_mesh
                # - NOT allowed by type mismatch => using just solps_extra._inside_mesh
                # => solps_extra_emission will be sampled where overlapping with original and where 0
                emission_2d_raw = solps_extra.halpha_total_radiation
                inside_mesh_original = solps_simulation._inside_mesh
                #cr = solps_extra.mesh.cr
                #cz = solps_extra.mesh.cz
                #for i in range(emission_2d_raw.shape[0]):
                #    for j in range(emission_2d_raw.shape[1]):
                 #       emission_2d_raw[i,j] *= (1.0 - inside_mesh_original(cr[i,j], cz[i,j]))
                solps_simulation = solps_extra
        elif use_homemade_emission is True:
            emission_2d_raw = solps_simulation["halpha_total_emission"].T

    elif cfg['plasma']['type_radiation'] == "halpha_mol_radiation":
        emission_2d_raw = solps_simulation.halpha_mol_radiation

    elif cfg['plasma']['type_radiation'] == "total_radiation":
        emission_2d_raw = solps_simulation.total_radiation

    emission_max = emission_2d_raw.max()

    L1st = emission_2d_raw / ComputeGradientNorm(cfg = cfg, b2fgmtry = b2fgmtry, SOLPSsim = solps_simulation, is_homemade = use_homemade_emission)
    L2nd = emission_2d_raw / ComputeHessianNorm2(cfg = cfg, b2fgmtry = b2fgmtry, SOLPSsim = solps_simulation, is_homemade = use_homemade_emission) / L1st

    # avoid pathetic behaviour of np.min and np.max when nan is present
    L1st[np.where(np.isnan(L1st))] = np.inf
    L2nd[np.where(np.isnan(L2nd))] = np.inf

    L1st /= alpha
    L2nd /= alpha

    fancyStep = np.zeros(L1st.shape)

    for iy in range(fancyStep.shape[0]):
        for ix in range(fancyStep.shape[1]):
            fancyStep[iy,ix] = np.max([step_min, np.min([step_max, L1st[iy,ix], L2nd[iy,ix]])])

    step_function_2d_raw = step_max * (fancyStep / step_max) ** (beta * (emission_2d_raw / emission_max))
    step_function_2d_raw[step_function_2d_raw < step_min] = step_min

    if use_homemade_emission is False:
        mesh = solps_simulation.mesh
        # if extra is False: inside_mesh = solps_simulation._inside_mesh
        inside_mesh = solps_simulation._inside_mesh

    #############################################################
    # padding to smooth transition from out- to in-side B2 mesh #
    #############################################################

    if cfg['raytracing']['sampling']['padding']['use_padding']:

        # load padded B2 mesh and over-write inside_mesh
        if extra is False:
            mesh = load_solps_mesh_padded(cfg = cfg)
        elif extra is True:
            mesh = load_solps_mesh_padded(cfg = cfg, extra = True)
        
        inside_mesh = SOLPSSimulation(mesh = mesh, species_list = ["i_am_useless"])._inside_mesh

        # load step_pad and compute step_pad_mean (log average)
        step_pad = cfg['raytracing']['sampling']['padding']['step_pad']
        step_pad_mean = pow(1E+01, np.mean([np.log10(step_max), np.log10(step_pad)]))
        steps_pad = [step_pad, step_pad_mean]

        # add padding to step_function_2d_raw and fill with step_pad*
        #
        # CAUTION. (ny, nx) CHERAB convention!
        
        for pad in steps_pad:
            padding_layer = np.ones((1, step_function_2d_raw.shape[1])) * pad
            step_function_2d_raw = np.concatenate([step_function_2d_raw, padding_layer], axis = 0)
            step_function_2d_raw = np.concatenate([padding_layer, step_function_2d_raw], axis = 0)

    #############################################################

    if cfg['raytracing']['sampling']['resize']['save'] is True:
        pickle.dump(step_function_2d_raw, open(os.path.join(cfg['input_directory'],
                    cfg['run_directory'], "step_function_2d_raw.pkl"), "wb"))

    if cfg['raytracing']['sampling']['resize']['load'] is True:
        step_function_2d_raw = pickle.load(open(os.path.join(cfg['input_directory'],
                               cfg['run_directory'], "step_function_2d_raw.pkl"), "rb"))

    ########################

    # create function

    # CAUTION.
    #
    # Outside B2 mesh (core & far-SOL>EIRENE) the emission is 0.0
    # => step = step_max in inhomogeneous.pyx
    #
    # Similarly for extra SOLPS emission => step = step_max where
    # overlapping with original SOLPS emission

    if use_homemade_emission is False:
        step_function_2d = SOLPSFunction2D.instance(inside_mesh, step_function_2d_raw)
        if extra is True: step_function_2d *= (1.0 - inside_mesh_original)
    elif use_homemade_emission is True:
        step_function_2d = Interpolate2DLinear(b2fgmtry['x_centres'], b2fgmtry['y_centres'], step_function_2d_raw.T, extrapolate = True)
    step_function_3d = AxisymmetricMapper(step_function_2d)

    ########################

    if plot_step_function is True:

        pad = 1E-03
        if use_homemade_emission is False:
            minr = mesh.mesh_extent['minr'] - pad
            maxr = mesh.mesh_extent['maxr'] + pad
            minz = mesh.mesh_extent['minz'] - pad
            maxz = mesh.mesh_extent['maxz'] + pad
        elif use_homemade_emission is True:
            minr = b2fgmtry['x_centres'][0] - pad
            maxr = b2fgmtry['x_centres'][-1] + pad
            minz = b2fgmtry['y_centres'][0] - pad
            maxz = b2fgmtry['y_centres'][-1] + pad

        plot_2d(ri = np.array([minr, maxr]), zi = np.array([minz, maxz]), f2d = step_function_2d, scale = "log", eq = None, title = r'Non-uniform step [m]')

        if cfg['plotting']['save_figures'] is True:
            plt.savefig(os.path.join(cfg['output_directory_extended'], 'step_function.png'))

    return step_function_3d