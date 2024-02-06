import os
import pickle
import numpy as np
from solps_python_scripts.clean import clean 
from solps_python_scripts.read_ft44_rfield import read_ft44_rfield
from solps_python_scripts.read_b2fgmtry import read_b2fgmtry
from solps_python_scripts.read_b2fstate import read_b2fstate
from solps_python_scripts.read_b2fplasmf import read_b2fplasmf
from solps_python_scripts.utilities.accessories import load_pickle

def read_ft44(where = ".", verbose = True, save = True):


    # Author: Wouter Dekeyser
    # November 2016
    #
    # Re-written in python by: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

    
     b2fstate  = read_b2fstate(where = where, verbose = False, save = True)
     b2fplasmf = read_b2fplasmf(where = where, verbose = False, save = True)
     neut = {}
     wld = {}   
     try:
        fid = open(os.path.join(where, "fort.44"), "r")

        neut = {}
        wld  = {}

        eV = 1.6022E-19 

        ## Read dimensions

        # nx, ny, version
        line = fid.readline().split()
        nx   = int(line[0])
        ny   = int(line[1])
        ver  = int(line[2])

        if ver != 20081111 and ver != 20160829 and ver != 20170328 and ver != 20180323 and ver != 20201006:
            raise ValueError('read_ft44: unknown format of fort.44 file')

        if verbose is True: print('read_ft44 -- file version ' + str(ver))

        # natm, nmol, nion
        line = fid.readline().split()
        natm = int(line[0])
        nmol = int(line[1])
        nion = int(line[2])
        nfla = int(b2fstate['ns'] - natm)
        
        neut['natm'] = natm
        neut['nmol'] = nmol
        neut['nion'] = nion
        neut['nfla'] = nfla

        neut['species'] = []

        # for now, ignore reading species labels
        atms = []
        mols = []
        ions = []
        for i in range(natm):
            atms += fid.readline().split()
        for i in range(nmol):
            mols += fid.readline().split()
        for i in range(nion):
            ions += fid.readline().split()

        neut['atms'] = atms
        neut['mols'] = mols
        neut['ions'] = ions

        #################################################################################


        # cumulative radiation emission
        
        neut['eneutrad'] = np.abs(read_ft44_rfield(fid,ver,'eneutrad',[nx,ny,natm]))
        print("trovato")
#        except: pass
        try: neut['emolrad']  = np.abs(read_ft44_rfield(fid,ver,'emolrad', [nx,ny,nmol]))
        except: pass
        try: neut['eionrad']  = np.abs(read_ft44_rfield(fid,ver,'eionrad', [nx,ny,nion]))
        except: pass
        try: neut['etotrad']  = neut['eneutrad'].sum(axis = 2) + neut['emolrad'].sum(axis = 2) + neut['eionrad'].sum(axis = 2)
        except: pass

        try: neut['totrad']   = neut['eneutrad'] + b2fplasmf['rqradtot'] + b2fplasmf['rqbrmtot']
        except: pass


        fid.close()

        if save is True:
            fort44 = {}
            fort44['neut'] = clean(neut)
            fort44['wld'] = clean(wld)
            pickle.dump(fort44, open(os.path.join(where, "fort.44.pkl"), "wb"))
     except: pass

     return (clean(neut), clean(wld))
