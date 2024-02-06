
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
import re
import pickle
import numpy as np
import scipy.io

from raysect.core import Point2D, Vector2D

from cherab.solps.mesh_geometry import SOLPSMesh
from cherab.solps import SOLPSSimulation

from solps_python_scripts.read_b2fgmtry  import read_b2fgmtry
from solps_python_scripts.read_b2fstate  import read_b2fstate
from solps_python_scripts.read_b2fplasmf import read_b2fplasmf
from solps_python_scripts.read_ft44      import read_ft44

def load_solps_data(cfg = None, extra = False):

    """
    Load a SOLPS simulation from SOLPS-ITER output files

    :param cfg.json
    :rtype: SOLPSSimulation

    #####################################################################

    MMM > TODO:

    - NOT all variables currently included
    - species_density:
        - b2fplasmf['na'] > FLUID NEUTRALS IN POSITION 0
        - *[:,:,0] = fort44['dab2'] > temporary work-around solution
    - sim.neutral_temperature =  fort44['tab2'] > put in eirene
    - add sim.halpha_atm_radiation field
    - sim.total_radiation > Halpha NOT included? Or already implicitly in?
    - sim.eirene:
        - update eirene.py to resemble B2.5 (*_f2d and *_f3d) IF POSSIBLE
        - insert neutral data there from fort.44 and fort.46
    - suggest to stick to SOLPS-ITER nomenclature
    - fix species_list: likely bug in CHERAB

    """

    ####################################
    # Load data from SOLPS-ITER output #
    ####################################

    SOLPS_directory = os.path.join(cfg['baserun'],
                                   'input',
                                   cfg['run'],
                                   cfg['plasma']['SOLPS']['SOLPS_directory'])

    if extra is True:
        SOLPS_directory = cfg['plasma']['SOLPS']['extra_directory']

    # try:

    #     b2fstate   = pickle.load(open(os.path.join(SOLPS_directory, 'b2fstate.pkl'), 'rb'))
    #     b2fplasmf   = pickle.load(open(os.path.join(SOLPS_directory, 'b2fplasmf.pkl'), 'rb'))
    #     fort44 = pickle.load(open(os.path.join(SOLPS_directory, 'fort44.pkl'), 'rb'))
    #     fort44 = fort44['neut']

    #     species = b2fstate['species']
    #     mesh    = load_solps_mesh(cfg = cfg)
    #     sim     = SOLPSSimulation(mesh = mesh, species_list = species)

    # except:

    try:
        b2fstate = read_b2fstate(where = SOLPS_directory, save = True)
    except:
        print()
        print('NO b2fstate => basic species assumed: ns = 2 ~ D0+ and D1+')
        print()
        b2fstate = {}
        b2fstate['ns'] = int(2)
        b2fstate['zamax'] = [0, 1]
        b2fstate['species'] = ["D0+", "D1+"]

    [neut, _] = read_ft44(where = SOLPS_directory, save = True)
    fort44 = neut

    # populate species (B2 convention):
    #
    # - all ionised stages
    # - fluid neutral atoms

    species = b2fstate['species']
    mesh    = load_solps_mesh(cfg = cfg, extra = extra)   
    sim     = SOLPSSimulation(mesh = mesh, species_list = species)

    # [_, b2fplasmf] = read_b2fplasmf(where = SOLPS_directory,
    #                                 nx = mesh.nx,
    #                                 ny = mesh.ny,
    #                                 ns = b2fstate['ns'],
    #                                 save = True)
    b2fplasmf = read_b2fplasmf(where = SOLPS_directory, save = True)

    # populate SOLSSimulation
    #
    # - SOLPS b2f*    convention: (nx+2) x (ny+2)
    # - SOLPS fort.4* convention:  nx    x  ny
    # - CHERAB        convention:  ny    x  nx

    # data must be C-continuous => BOH

    # temperatures [eV]

    te = b2fplasmf['te'].T.copy(order = 'C')
    ti = b2fplasmf['ti'].T.copy(order = 'C')

    sim.electron_temperature = te
    sim.ion_temperature =      ti
    # sim.neutral_temperature =  fort44['tab2'].T

    # densities [m^{-3}]
#    
#    ne = b2fplasmf['ne'].T.copy(order = 'C')
#    na = b2fplasmf['na'].T.copy(order = 'C')
#    dab2 = fort44['dab2'].T.copy(order = 'C')
#    j = 0 # populate neutral ATOMS on B2 mesh
#    for i in range(b2fstate['ns']):
#        if b2fstate['zamax'][i] == 0:
#            na[i,:,:] = dab2[j,:,:]
#            j += 1
#
#    sim.electron_density = ne
#    sim.species_density =  na
#
#    # Halpha radiation [ph / m^3 / s]
#
#    emissmol = np.abs(fort44['emissmol'].sum(axis = 2).T.copy(order = 'C'))
#    emiss    = np.abs(fort44['emiss'].sum(axis = 2).T.copy(order = 'C'))
#
#    sim.halpha_mol_radiation =   emissmol
#    sim.halpha_total_radiation = emissmol + emiss

    ##########################################################################

    # MMM: temporary ad hoc fixes

    # MMM: different weights

    use_extra_emission    = cfg['plasma']['SOLPS']['use_extra_emission']
    use_weighted_emission = cfg['plasma']['SOLPS']['weights']['use_weighted_emission']
    if use_weighted_emission is True:
        if use_extra_emission is False: raise ValueError('Extra emission needed to use weights!')
        if   extra is False: weight = cfg['plasma']['SOLPS']['weights']['weight_main_emission']
        elif extra is True:  weight = cfg['plasma']['SOLPS']['weights']['weight_extra_emission']

    # MMM: fix_artifact
    if cfg['plasma']['SOLPS']['artifact']['fix_artifact'] is True:

        start_index = cfg['plasma']['SOLPS']['artifact']['start_index']
        end_index = cfg['plasma']['SOLPS']['artifact']['end_index']
        threshold_value = cfg['plasma']['SOLPS']['artifact']['threshold_value']
        new_value = cfg['plasma']['SOLPS']['artifact']['new_value']

        tmp = emissmol + emiss
        for iy in range(tmp.shape[0]):
            for ix in range(tmp.shape[1]):
                if ix > start_index and end_index < 60:
                    if tmp[iy,ix] > threshold_value:
                        tmp[iy,ix] = new_value
        sim.halpha_total_radiation = tmp

    ##########################################################################

    # CAUTION.
    #
    # - Bremsstrahlung + Line radiation from {Z>1-ions ; atoms ; molecules ; molecular ions}
    # - rq* must be summed all over species (although, if EIRENE = ON, fluid neutrals do NOT contribute)
    # - rq*   now in units [W] coherently with SOLPS-ITER convention
    # - e*rad now in units [W] coherently with SOLPS-ITER convention

    try:    rqbrm = np.abs(b2fplasmf['rqbrm'].sum(axis = 2).T.copy(order = 'C'))
    except: rqbrm = np.abs(b2fplasmf['rqbrm'].T.copy(order = 'C'))
    try:    rqrad = np.abs(b2fplasmf['rqrad'].sum(axis = 2).T.copy(order = 'C'))
    except: rqrad = np.abs(b2fplasmf['rqrad'].T.copy(order = 'C'))

    try: eneutrad = np.abs(fort44['eneutrad'].sum(axis = 2).T.copy(order = 'C'))
    except: pass
    try: emolrad  = np.abs(fort44['emolrad'].sum(axis = 2).T.copy(order = 'C'))
    except: pass
    try: eionrad  = np.abs(fort44['eionrad'].sum(axis = 2).T.copy(order = 'C'))
    except: pass
    
    if cfg['plasma']['SOLPS']['rqbrm'] is True:
       sim.total_radiation  = rqbrm         # [W]
       if  cfg['plasma']['SOLPS']['rqrad'] is True:
           sim.total_radiation  += rqrad        # [W]

    if cfg['plasma']['SOLPS']['rqbrm'] is False:
        sim.total_radiation  = rqrad         # [W]

    try: sim.total_radiation += eneutrad      # [W]
    except: pass
    try: sim.total_radiation += emolrad  # [W]
    except: pass
    try: sim.total_radiation += eionrad  # [W]
    except: pass
    
    sim.total_radiation /= mesh.vol # [W / m^3]
    
    if use_weighted_emission is True:
        sim.total_radiation        *= weight
        sim.halpha_total_radiation *= weight

    sim.eirene = None

    return sim

def load_solps_mesh(cfg = None, extra = False):

    """
    Load the SOLPS mesh geometry from SOLPS-ITER output files
    """

    SOLPS_directory = os.path.join(cfg['baserun'],
                                   'input',
                                   cfg['run'],
                                   cfg['plasma']['SOLPS']['SOLPS_directory'])

    if extra is True:
        SOLPS_directory = cfg['plasma']['SOLPS']['extra_directory']

    # try: gmtry = pickle.load(open(os.path.join(SOLPS_directory, 'b2fgmtry.pkl'), 'rb'))
    # except: gmtry = read_b2fgmtry(where = SOLPS_directory, save = True)
    gmtry = read_b2fgmtry(where = SOLPS_directory, save = True)
    
    vr = gmtry['crx']
    vz = gmtry['cry']
    vol = gmtry['vol']

    lix = gmtry['leftix']
    rix = gmtry['rightix']
    tix = gmtry['topix']
    bix = gmtry['bottomix']

    liy = gmtry['leftiy']
    riy = gmtry['rightiy']
    tiy = gmtry['topiy']
    biy = gmtry['bottomiy']

    # MMM: speed-up factor scaling
    if cfg['raytracing']['sampling']['resize']['use_resize'] is True:
          rescale = cfg['plasma']['sampling']['resize']['resize_factor']
    else: rescale = 1.0

    vz = vz * rescale
    vr = vr * rescale
    vol = vol * rescale**3

    lix = lix
    rix = rix
    tix = tix
    bix = bix

    liy = liy
    riy = riy
    tiy = tiy
    biy = biy

    # building ix neighbors: watch out the indeces
    neighbix = np.zeros((4, (np.shape(lix))[1], (np.shape(lix))[0]))
    neighbix[0,:,:] = lix.T
    neighbix[1,:,:] = bix.T
    neighbix[2,:,:] = rix.T
    neighbix[3,:,:] = tix.T

    # building iy neighbors: watch out the indeces
    neighbiy = np.zeros((4, (np.shape(liy))[1], (np.shape(liy))[0]))
    neighbiy[0,:,:] = liy.T
    neighbiy[1,:,:] = biy.T
    neighbiy[2,:,:] = riy.T
    neighbiy[3,:,:] = tiy.T

    # build mesh object
    mesh = SOLPSMesh(vr.T,
                     vz.T,
                     vol.T,
                     neighbix,
                     neighbiy)

    # Load the ***CENTRE*** points of the grid cells.
    cr = mesh._cr.T
    cz = mesh._cz.T

    # Load cell basis vectors
    nx = mesh.nx
    ny = mesh.ny

    # 3-index array/list/boh => need to understand what 2 stands for
    cell_poloidal_basis = np.empty((nx, ny, 2), dtype = object)

    # now you use the CENTRE points of the RECTANGULAR cells (cr and cz)
    # to do strange things....
    #
    for i in range(nx):
        for j in range(ny):

            # Work out cell's 2D parallel vector in the poloidal plane:
            if i == nx - 1:
                # Special case for end of array, repeat previous calculation.
                # This is because I don't have access to the gaurd cells.
                xp_x = cr[i, j] - cr[i-1, j]
                xp_y = cz[i, j] - cz[i-1, j]
                norm = np.sqrt(xp_x**2 + xp_y**2)
                cell_poloidal_basis[i, j, 0] = Point2D(xp_x/norm, xp_y/norm)
            else:
                xp_x = cr[i+1, j] - cr[i, j]
                xp_y = cz[i+1, j] - cz[i, j]
                norm = np.sqrt(xp_x**2 + xp_y**2)
                cell_poloidal_basis[i, j, 0] = Point2D(xp_x/norm, xp_y/norm)

            # Work out cell's 2D radial vector in the poloidal plane
            if j == ny - 1:
                # Special case for end of array, repeat previous calculation.
                yr_x = cr[i, j] - cr[i, j-1]
                yr_y = cz[i, j] - cz[i, j-1]
                norm = np.sqrt(yr_x**2 + yr_y**2)
                cell_poloidal_basis[i, j, 1] = Point2D(yr_x/norm, yr_y/norm)
            else:
                yr_x = cr[i, j+1] - cr[i, j]
                yr_y = cz[i, j+1] - cz[i, j]
                norm = np.sqrt(yr_x**2 + yr_y**2)
                cell_poloidal_basis[i, j, 1] = Point2D(yr_x/norm, yr_y/norm)

    mesh._poloidal_grid_basis = cell_poloidal_basis

    return mesh

#########################################################################################
#########################################################################################
#########################################################################################

def load_solps_mesh_padded(cfg = None, extra = False):

    """
    Load the SOLPS mesh geometry from SOLPS-ITER output files but add padding on N and S
    """

    SOLPS_directory = os.path.join(cfg['baserun'],
                                   'input',
                                   cfg['run'],
                                   cfg['plasma']['SOLPS']['SOLPS_directory'])

    if extra is True:
        SOLPS_directory = cfg['plasma']['SOLPS']['extra_directory']

    # try: gmtry = pickle.load(open(os.path.join(SOLPS_directory, 'b2fgmtry.pkl'), 'rb'))
    # except: gmtry = read_b2fgmtry(where = SOLPS_directory, save = True)
    gmtry = read_b2fgmtry(where = SOLPS_directory, save = True)
    
    vr = gmtry['crx']
    vz = gmtry['cry']
    vol = gmtry['vol']

    lix = gmtry['leftix']
    rix = gmtry['rightix']
    tix = gmtry['topix']
    bix = gmtry['bottomix']

    liy = gmtry['leftiy']
    riy = gmtry['rightiy']
    tiy = gmtry['topiy']
    biy = gmtry['bottomiy']

    # MMM: speed-up factor scaling
    if cfg['raytracing']['sampling']['resize']['use_resize'] is True:
        rescale = cfg['raytracing']['sampling']['resize']['resize_factor']
    else: rescale = 1E+00

    vz = vz * rescale
    vr = vr * rescale
    vol = vol * rescale**3

    lix = lix
    rix = rix
    tix = tix
    bix = bix

    liy = liy
    riy = riy
    tiy = tiy
    biy = biy

    #####################################################################

    ###########
    # PADDING #
    ###########

    # - B2 mesh with NO ghost cells at this stage
    # - New cell layers (x2) added radially to existing B2 mesh
    # - Radially-pointing vector aligned with y-sides of the cells
    # - Volume/neighbor data not needed <=> NO emission in padding layers
    #
    # - Sampling step in padding layers:
    #   - Outermost = step_pad_mean = log average between step_max and step_pad
    #   - Innermost = step_pad
    # => step_pad in innermost determine the resolution of the discontinuity
    #
    # - Width of padding layers:
    #   - Outermost = step_max => at least one sample point with step_max
    #   - Innermost = step_pad => at least one sample point with step_pad_mean
    #
    # CAUTION.
    #
    # Not as good when sampling EXITS B2 mesh => discontinuity not necessarily well resolved
    # but anyway a bit ameliorated!!

    nx       = vr.shape[0] 
    npad     = 2
    step_pad = cfg['raytracing']['sampling']['padding']['step_pad']
    step_max = cfg['raytracing']['sampling']['step_max']

    # log average
    step_pad_mean = pow(1E+01, np.mean([np.log10(step_max), np.log10(step_pad)]))

    # - index 0: width = step_pad_mean, i.e. first  layer outside B2 mesh
    # - index 1: width = step_max,      i.e. second layer outside B2 mesh
    width_pad = [step_pad_mean, step_max]

    newSr = np.zeros((nx, npad, 4))
    newSz = np.zeros((nx, npad, 4))
    newNr = np.zeros((nx, npad, 4))
    newNz = np.zeros((nx, npad, 4))

    dummy = np.zeros((nx, 1))

    for ipad in range(npad):

        for ix in range(nx):

            # CAUTION.
            # - N and S boundaries only
            # - If ipad == 0: iy = 0 and iy = vr.shape[1]-1 are B2 original cells
            # - If ipad >= 1: iy = 0 and iy = vr.shape[1]-1 are the new padding layers

            for iy in [0, vr.shape[1]-1]:

                pts = []
                for iv in range(4): pts += [Point2D(vr[ix,iy,iv], vz[ix,iy,iv])]

                # CAUTION.
                # - Assumed +y direction
                # - B2 vertex ordering convention
                
                v02 = pts[0].vector_to(pts[2]).normalise()
                v13 = pts[1].vector_to(pts[3]).normalise()

                # S boundary
                if iy == 0:

                    # same vertices as on S boundary
                    newSr[ix,ipad,2] = vr[ix,iy,0]
                    newSz[ix,ipad,2] = vz[ix,iy,0]
                    newSr[ix,ipad,3] = vr[ix,iy,1]
                    newSz[ix,ipad,3] = vz[ix,iy,1]

                    # inward pointing
                    v02 *= -1
                    v13 *= -1

                    # new vertices
                    newSr[ix,ipad,0] = (pts[0] + v02 * width_pad[ipad]).x
                    newSz[ix,ipad,0] = (pts[0] + v02 * width_pad[ipad]).y
                    newSr[ix,ipad,1] = (pts[1] + v13 * width_pad[ipad]).x
                    newSz[ix,ipad,1] = (pts[1] + v13 * width_pad[ipad]).y                      

                # N boundary
                elif iy == vr.shape[1]-1:

                    # same vertices as on N boundary
                    newNr[ix,ipad,0] = vr[ix,iy,2]
                    newNz[ix,ipad,0] = vz[ix,iy,2]
                    newNr[ix,ipad,1] = vr[ix,iy,3]
                    newNz[ix,ipad,1] = vz[ix,iy,3]

                    # already outward pointing

                    # new vertices
                    newNr[ix,ipad,2] = (pts[2] + v02 * width_pad[ipad]).x
                    newNz[ix,ipad,2] = (pts[2] + v02 * width_pad[ipad]).y
                    newNr[ix,ipad,3] = (pts[3] + v13 * width_pad[ipad]).x
                    newNz[ix,ipad,3] = (pts[3] + v13 * width_pad[ipad]).y  

        # S boundary                

        # add new TOP (index = 0 ~ S) layer to vr and vz
        padding_layer = np.reshape(newSr[:,ipad,:], (nx, 1, 4))
        vr = np.concatenate([padding_layer, vr], axis = 1)
        padding_layer = np.reshape(newSz[:,ipad,:], (nx, 1, 4))
        vz = np.concatenate([padding_layer, vz], axis = 1)

        # add dummy TOP (index = 0 ~ S) layer to vr and vz
        vol = np.concatenate([dummy, vol], axis = 1)
        lix = np.concatenate([dummy, lix], axis = 1)
        rix = np.concatenate([dummy, rix], axis = 1)
        tix = np.concatenate([dummy, tix], axis = 1)
        bix = np.concatenate([dummy, bix], axis = 1)
        liy = np.concatenate([dummy, liy], axis = 1)
        riy = np.concatenate([dummy, riy], axis = 1)
        tiy = np.concatenate([dummy, tiy], axis = 1)
        biy = np.concatenate([dummy, biy], axis = 1)

        # N boundary

        # add new BOTTOM (index = ny ~ N) layer to vr and vz
        padding_layer = np.reshape(newNr[:,ipad,:], (nx, 1, 4))
        vr = np.concatenate([vr, padding_layer], axis = 1)
        padding_layer = np.reshape(newNz[:,ipad,:], (nx, 1, 4))
        vz = np.concatenate([vz, padding_layer], axis = 1)

        # add dummy BOTTOM (index = ny ~ N) layer to vr and vz
        vol = np.concatenate([vol, dummy], axis = 1)
        lix = np.concatenate([lix, dummy], axis = 1)
        rix = np.concatenate([rix, dummy], axis = 1)
        tix = np.concatenate([tix, dummy], axis = 1)
        bix = np.concatenate([bix, dummy], axis = 1)
        liy = np.concatenate([liy, dummy], axis = 1)
        riy = np.concatenate([riy, dummy], axis = 1)
        tiy = np.concatenate([tiy, dummy], axis = 1)
        biy = np.concatenate([biy, dummy], axis = 1)

    #####################################################################

    # building ix neighbors: watch out the indeces
    neighbix = np.zeros((4, (np.shape(lix))[1], (np.shape(lix))[0]))
    neighbix[0,:,:] = lix.T
    neighbix[1,:,:] = bix.T
    neighbix[2,:,:] = rix.T
    neighbix[3,:,:] = tix.T

    # building iy neighbors: watch out the indeces
    neighbiy = np.zeros((4, (np.shape(liy))[1], (np.shape(liy))[0]))
    neighbiy[0,:,:] = liy.T
    neighbiy[1,:,:] = biy.T
    neighbiy[2,:,:] = riy.T
    neighbiy[3,:,:] = tiy.T

    # build mesh object
    mesh = SOLPSMesh(vr.T,
                     vz.T,
                     vol.T,
                     neighbix,
                     neighbiy)

    # Load the ***CENTRE*** points of the grid cells.
    cr = mesh._cr.T
    cz = mesh._cz.T

    # Load cell basis vectors
    nx = mesh.nx
    ny = mesh.ny

    # 3-index array/list/boh => need to understand what 2 stands for
    cell_poloidal_basis = np.empty((nx, ny, 2), dtype = object)

    # now you use the CENTRE points of the RECTANGULAR cells (cr and cz)
    # to do strange things....
    #
    for i in range(nx):
        for j in range(ny):

            # Work out cell's 2D parallel vector in the poloidal plane:
            if i == nx - 1:
                # Special case for end of array, repeat previous calculation.
                # This is because I don't have access to the gaurd cells.
                xp_x = cr[i, j] - cr[i-1, j]
                xp_y = cz[i, j] - cz[i-1, j]
                norm = np.sqrt(xp_x**2 + xp_y**2)
                cell_poloidal_basis[i, j, 0] = Point2D(xp_x/norm, xp_y/norm)
            else:
                xp_x = cr[i+1, j] - cr[i, j]
                xp_y = cz[i+1, j] - cz[i, j]
                norm = np.sqrt(xp_x**2 + xp_y**2)
                cell_poloidal_basis[i, j, 0] = Point2D(xp_x/norm, xp_y/norm)

            # Work out cell's 2D radial vector in the poloidal plane
            if j == ny - 1:
                # Special case for end of array, repeat previous calculation.
                yr_x = cr[i, j] - cr[i, j-1]
                yr_y = cz[i, j] - cz[i, j-1]
                norm = np.sqrt(yr_x**2 + yr_y**2)
                cell_poloidal_basis[i, j, 1] = Point2D(yr_x/norm, yr_y/norm)
            else:
                yr_x = cr[i, j+1] - cr[i, j]
                yr_y = cz[i, j+1] - cz[i, j]
                norm = np.sqrt(yr_x**2 + yr_y**2)
                cell_poloidal_basis[i, j, 1] = Point2D(yr_x/norm, yr_y/norm)

    mesh._poloidal_grid_basis = cell_poloidal_basis

    return mesh
