
# Copyright 2016-2018 Euratom
# Copyright 2016-2018 United Kingdom Atomic Energy Authority
# Copyright 2016-2018 Centro de Investigaciones Energéticas, Medioambientales y Tecnológicas
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

import os
import pickle
import numpy as np

from raysect.core import Point2D, Vector2D

from eirene_simulation import EIRENESimulation

from solps_python_scripts.read_triangle_mesh import read_triangle_mesh
from solps_python_scripts.read_ft46          import read_ft46

from quantum_data import FractionExcitedHydrogenGas

def load_eirene_data(cfg = None):

    """
    Load a EIRENE simulation from SOLPS-ITER>EIRENE output files

    :param cfg.json
    :rtype: EIRENESimulation

    #####################################################################

    MMM > TODO:

    - put in SOLPSSimulation.eirene
    - add velocities (see read_ft46.py)

    """

    ################################
    # Load data from EIRENE output #
    ################################

    SOLPS_directory = os.path.join(cfg['baserun'],
                                   'input',
                                   cfg['run'],
                                   cfg['plasma']['SOLPS']['SOLPS_directory'])

    # MMM: different weights
    use_weighted_emission = cfg['plasma']['SOLPS']['weights']['use_weighted_emission']
    if use_weighted_emission is True:
      # only weight of main emission by assumption (EIRENE from main emission)
      weight = cfg['plasma']['SOLPS']['weights']['weight_main_emission']

    try:
      mesh   = pickle.load(open(os.path.join(SOLPS_directory, 'triangle_mesh.pkl'), 'rb'))
    except:
      mesh   = read_triangle_mesh(where = SOLPS_directory, save = True)
    fort46 = read_ft46(where = SOLPS_directory, save = True)
  
    sim = EIRENESimulation(mesh = mesh)

    # populate EIRENESimulation

    # densities [m^{-3}]
    
#    sim.pdena = fort46['pdena'][:,0] # only main ion for Halpha
#    sim.pdenm = fort46['pdenm']
#    sim.pdeni = fort46['pdeni']
#
#    # energy densities [Pa]
#    
#    sim.edena = fort46['edena'][:,0] # only main ion for Halpha
#    sim.edenm = fort46['edenm']
#    sim.edeni = fort46['edeni']
#
#    # temperatures [eV]
#
#    sim.tdena = fort46['tdena'][:,0] # only main ion for Halpha
#    sim.tdenm = fort46['tdenm']
#    sim.tdeni = fort46['tdeni']

    # Ha emission [photons / m^{3} / s]
    #
    # CAUTION.
    #
    # - ref: [https://iopscience.iop.org/article/10.1088/1361-6587/abffb7/meta#fnref-ppcfabffb7bib15]
    # - NO molecules nor molecular ions???

    # Boltzmann-like (i.e. neutral gas) fraction of excited hydrogen
    f_3  = FractionExcitedHydrogenGas(n_quantum = 3,
                                      T_eV = sim.tdena)

    sim.Ha_emission = sim.Ha_Einstein_coeff * f_3 * sim.pdena

    if use_weighted_emission is True: sim.Ha_emission *= weight

    return sim