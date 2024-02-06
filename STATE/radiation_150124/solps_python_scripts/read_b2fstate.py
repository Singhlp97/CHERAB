import os
import pickle
import numpy as np
import mendeleev as mdv
from solps_python_scripts.clean import clean
from solps_python_scripts.delete_ghost_cells import delete_ghost_cells
from solps_python_scripts.read_rfield import read_rfield
from solps_python_scripts.read_ifield import read_ifield
from solps_python_scripts.utilities.accessories import load_pickle

def read_species(b2fstate = None):

    species = []

    elements = mdv.get_all_elements()

    for i in range(b2fstate['ns']):
        am = int(b2fstate['am'][i])
        zn = int(b2fstate['zn'][i])
        za = int(b2fstate['zamax'][i])
        for element in elements:
            if element.atomic_number == zn:
                # deuterium (for coherence with fort.44 notation)
                if zn == 1 and am == 2:
                    species += ['D' + str(za) + "+"]
                #tritium (for coherence with fort.44 notation)
                elif zn == 1 and am == 3:
                    species += ['T' + str(za) + "+"]
                # other element
                else:
                    species += [element.symbol + str(za) + "+"]

    return species

def read_b2fstate(where = ".", verbose = True, save = True):

    # state = read_b2fstate(file)
    #
    # Read b2fstati/b2fstate file created by B2.5.
    #
    # Output is a struct "state" with all the data fields in the b2fstate/i
    # file (except for the dimensions nx,ny,ns, which are implicit in the array
    # sizes).
    #

    # Author: Wouter Dekeyser
    # November 2016
    #
    # Re-written in python by: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

    try:        
        
        state = load_pickle(where = where, what = "b2fstate", verbose = verbose)

    except:

        # Open file
        try:
            fid = open(os.path.join(where, "b2fstate"), "r")
        except:
            fid = open(os.path.join(where, "b2fstati"), "r")
            print()
            print('read_b2fstate -- b2fstati read instead')
            print()


        ## Get version of the b2fstate file

        line    = fid.readline()
        version = line[7:17]

        state = {}

        eV = 1.6022E-19

        if verbose is True: print('read_b2fstate -- file version ' + version)
        state['version'] = version


        ## Read dimensions nx, ny, ns

        dim = read_ifield(fid,'nx,ny,ns',3)
        nx  = int(dim[0])
        ny  = int(dim[1])
        ns  = int(dim[2])
        
        state['ns'] = ns

        fluxdim  = [nx+2,ny+2,2]
        fluxdimp = [nx+2,ny+2]
        fluxdims = [nx+2,ny+2,2,ns]
        if version >= '03.001.000':
            fluxdim  = [nx+2,ny+2,2,2]
            fluxdimp = fluxdim
            fluxdims = [nx+2,ny+2,2,2,ns]

        ## Read charges etc.

        state['zamin'] = read_rfield(fid,'zamin',ns)
        state['zamax'] = read_rfield(fid,'zamax',ns)
        state['zn']    = read_rfield(fid,'zn   ',ns)
        state['am']    = read_rfield(fid,'am   ',ns)

        state['species'] = read_species(b2fstate = state)

        ## Read state variables

        state['na']     = read_rfield(fid,'na'    ,[nx+2,ny+2,ns])
        state['ne']     = read_rfield(fid,'ne'    ,[nx+2,ny+2])
        state['ua']     = read_rfield(fid,'ua'    ,[nx+2,ny+2,ns])
        state['uadia']  = read_rfield(fid,'uadia' ,[nx+2,ny+2,2,ns])
        state['te']     = read_rfield(fid,'te'    ,[nx+2,ny+2]) / eV
        state['ti']     = read_rfield(fid,'ti'    ,[nx+2,ny+2]) / eV
        state['po']     = read_rfield(fid,'po'    ,[nx+2,ny+2])

        ## add nae and Zeff

        state['nae']  = np.zeros(state['na'].shape)
        state['Zeff'] = np.zeros(state['ne'].shape)
        for isp in range(ns):
            state['nae'][:,:,isp] = state['na'][:,:,isp] / state['ne']
            state['Zeff'] += state['zamin'][isp]**2 * state['nae'][:,:,isp]



        ## Close file

        fid.close()

        state = delete_ghost_cells(state)

        if save is True:
            pickle.dump(clean(state), open(os.path.join(where, "b2fstate.pkl"), "wb"))

    return clean(state)
