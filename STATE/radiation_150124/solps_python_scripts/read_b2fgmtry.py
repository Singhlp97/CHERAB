import os
import pickle
import numpy as np
from solps_python_scripts.clean import clean
from solps_python_scripts.delete_ghost_cells import delete_ghost_cells
from solps_python_scripts.read_ifield import read_ifield
from solps_python_scripts.read_rfield import read_rfield
from solps_python_scripts.utilities.accessories import load_pickle

def read_b2fgmtry(where = ".", verbose = True, save = True):

    # gmtry = read_b2fgmtry(file)
    #
    # Read b2fgmtry file created by B2.5.
    #
    #

    # Author: Wouter Dekeyser
    # November 2016
    #
    # Re-written in python by: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

    # CAUTION.
    # - Search b2fgmtry where indicated by the user
    # - If NO b2fgmtry there, search in ./baserun
    # - If NO b2fgmtry there, search in ../baserun
    # - If NO b2fgmtry there, error 

    try:

        gmtry = load_pickle(where = where, what = "b2fgmtry", verbose = verbose)

    except:

        if os.path.exists(os.path.join(where, "b2fgmtry")) is True:
            fid = open(os.path.join(where, "b2fgmtry"), "r")
        elif os.path.exists(os.path.join(where, "../baserun/b2fgmtry")) is True:
            fid = open(os.path.join(where, "../baserun/b2fgmtry"), "r")
        else: raise ValueError('b2fgmtry can NOT be found :(')

        ## Get version of the b2fstate file

        line    = fid.readline()
        version = line[7:17]

        if verbose is True: print("read_b2fgmtry -- file version " + version)

        ## Read dimensions nx, ny, and symmetry

        dim = read_ifield(fid,'nx,ny',[2])
        nx  = int(dim[0])
        ny  = int(dim[1])

        # Expected array sizes, gmtry
        qcdim = [nx+2, ny+2]
        if version >= '03.001.000':
            qcdim = [nx+2, ny+2, 2]

        gmtry = {}
        
        gmtry['nx'] = nx
        gmtry['ny'] = ny

        ## Read symmetry information

        gmtry['isymm'] = read_ifield(fid,'isymm',1)


        ## Read gmtry variables

        gmtry['crx']  = read_rfield(fid,'crx' ,[nx+2,ny+2,4])
        gmtry['cry']  = read_rfield(fid,'cry' ,[nx+2,ny+2,4])
        gmtry['vol']  = read_rfield(fid,'vol' ,[nx+2,ny+2])
        
        gmtry['nncut']     = read_ifield(fid,'nncut'    ,1)
        gmtry['leftcut']   = read_ifield(fid,'leftcut'  ,gmtry['nncut'])
        gmtry['rightcut']  = read_ifield(fid,'rightcut' ,gmtry['nncut'])
        gmtry['topcut']    = read_ifield(fid,'topcut'   ,gmtry['nncut'])
        gmtry['bottomcut'] = read_ifield(fid,'bottomcut',gmtry['nncut'])

        gmtry['leftix']    = read_ifield(fid,'leftix'   ,[nx+2,ny+2])
        gmtry['rightix']   = read_ifield(fid,'rightix'  ,[nx+2,ny+2])
        gmtry['topix']     = read_ifield(fid,'topix'    ,[nx+2,ny+2])
        gmtry['bottomix']  = read_ifield(fid,'bottomix' ,[nx+2,ny+2])
        gmtry['leftiy']    = read_ifield(fid,'leftiy'   ,[nx+2,ny+2])
        gmtry['rightiy']   = read_ifield(fid,'rightiy'  ,[nx+2,ny+2])
        gmtry['topiy']     = read_ifield(fid,'topiy'    ,[nx+2,ny+2])
        gmtry['bottomiy']  = read_ifield(fid,'bottomiy' ,[nx+2,ny+2])

        ## Close file

        fid.close()

        gmtry = delete_ghost_cells(gmtry)

        if save is True:
            pickle.dump(clean(gmtry), open(os.path.join(where, "b2fgmtry.pkl"), "wb"))

    return clean(gmtry)
