import os
import pickle
import numpy as np
from solps_python_scripts.clean import clean 
from solps_python_scripts.delete_ghost_cells import delete_ghost_cells
from solps_python_scripts.read_ifield import read_ifield
from solps_python_scripts.read_rfield import read_rfield
from solps_python_scripts.read_b2fgmtry import read_b2fgmtry
from solps_python_scripts.read_b2fstate import read_b2fstate
from solps_python_scripts.utilities.accessories import load_pickle

def read_b2fplasmf(where = ".", verbose = True, save = None):

    # [gmtry,state] = read_b2fplasmf(file,nx,ny,ns)
    #
    # Read formatted b2fplasmf file created by B2.5.
    #
    #

    # Author: Wouter Dekeyser
    # November 2016
    #
    # Re-written in python by: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

    try:

        state = load_pickle(where = where, what = "b2fplasmf")

    except:

        fid = open(os.path.join(where, "b2fplasmf"), "r")

        b2fgmtry = read_b2fgmtry(where = where, verbose = False)
        b2fstate = read_b2fstate(where = where, verbose = False)
        
        nx = b2fgmtry['nx']
        ny = b2fgmtry['ny']
        ns = b2fstate['ns']

        # gmtry can be loaded via read_b2fgmtry()

        # gmtry = {}
        state = {}

        eV = 1.6022E-19

        ## Get version of the b2fstate file

        line    = fid.readline()
        version = line[7:17]

        if verbose is True: print('read_b2fplasmf -- file version ' + version)

        # Expected array sizes, gmtry
        qcdim = [nx+2,ny+2]
        if version >= '03.001.000':
            qcdim  = [nx+2,ny+2,2]


        # Expected array sizes, state
        fluxdim  = [nx+2,ny+2,2]
        fluxdims = [nx+2,ny+2,2,ns]
        if version >= '03.001.000':
            fluxdim  = [nx+2,ny+2,2,2]
            fluxdims = [nx+2,ny+2,2,2,ns]

        ## Read gmtry variables



        ## Read state variables

        state['fna']    = read_rfield(fid,'fna'   ,fluxdims)
        state['fne']    = read_rfield(fid,'fne'   ,fluxdim)
        state['fni']    = read_rfield(fid,'fni'   ,fluxdim)
        state['na']     = read_rfield(fid,'na'    ,[nx+2,ny+2,ns])
        state['ne']     = read_rfield(fid,'ne'    ,[nx+2,ny+2])
        state['ne0']    = read_rfield(fid,'ne0'   ,[nx+2,ny+2])
        state['ne2']    = read_rfield(fid,'ne2'   ,[nx+2,ny+2])
        state['nep']    = read_rfield(fid,'nep'   ,[nx+2,ny+2])

        ## add nae and Zeff

        state['nae']  = np.zeros(state['na'].shape)
        state['Zeff'] = np.zeros(state['ne'].shape)
        for isp in range(ns):
            state['nae'][:,:,isp] = state['na'][:,:,isp] / state['ne']
            state['Zeff'] += b2fstate['zamin'][isp]**2 * state['nae'][:,:,isp]
        
        
        state['te']     = read_rfield(fid,'te'    ,[nx+2,ny+2]) / eV
        state['ti']     = read_rfield(fid,'ti'    ,[nx+2,ny+2]) / eV
        state['rqrad']        = read_rfield(fid,'rqrad'       ,[nx+2,ny+2,ns])
        state['rqbrm']        = read_rfield(fid,'rqbrm'       ,[nx+2,ny+2,ns])

        state['rqradtot'] = state['rqrad'].sum(axis = 2)
        state['rqbrmtot'] = state['rqbrm'].sum(axis = 2)


        state = delete_ghost_cells(state)

        if save is True:
            pickle.dump(clean(state), open(os.path.join(where, "b2fplasmf.pkl"), "wb"))

    return clean(state)
