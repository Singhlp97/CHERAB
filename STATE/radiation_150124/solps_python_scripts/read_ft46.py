import os
import pickle
import numpy as np
from solps_python_scripts.clean import clean 
from solps_python_scripts.read_ft44_rfield import read_ft44_rfield 
from solps_python_scripts.utilities.accessories import load_pickle

def read_ft46(where = ".", verbose = True, save = True):

    # tdata = read_ft46(file)
    #
    # Read fort.46 file. Convert to SI units.
    #
    # For now, only fort.46 version 20160513 recognized 
    #

    # Author: Wouter Dekeyser
    # November 2016
    #
    # Re-written in python by: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

    try:

        tdata = load_pickle(where = where, what = "fort.46")

    except:
        
        fid = open(os.path.join(where, "fort.46"), "r")

        ## Read dimensions

        # ntri, version
        line = fid.readline().split()
        ntri = int(line[0])
        ver  = int(line[1])

        if ver != 20160513 and ver != 20160829 and ver != 20170930:
            raise ValueError('read_ft46: unknown format of fort.46 file')

        if verbose is True: print("read_ft46: -- file version " + str(ver))

        # natm, nmol, nion
        line = fid.readline().split()
        natm = int(line[0])
        nmol = int(line[1])
        nion = int(line[2])

        # for now, ignore reading species labels
        for i in range(natm):
            line = fid.readline()
        for i in range(nmol):
            line = fid.readline()
        for i in range(nion):
            line = fid.readline()

        eV = 1.6022E-19

        ## Read data

        tdata = {}

        tdata['pdena']  = read_ft44_rfield(fid,ver,'pdena',[ntri,natm])*1e6 # m^{-3}
        tdata['pdenm']  = read_ft44_rfield(fid,ver,'pdenm',[ntri,nmol])*1e6
        tdata['pdeni']  = read_ft44_rfield(fid,ver,'pdeni',[ntri,nion])*1e6

        tdata['edena']  = read_ft44_rfield(fid,ver,'edena',[ntri,natm])*1e6*eV # J m^{-3}
        tdata['edenm']  = read_ft44_rfield(fid,ver,'edenm',[ntri,nmol])*1e6*eV
        tdata['edeni']  = read_ft44_rfield(fid,ver,'edeni',[ntri,nion])*1e6*eV

        tdata['tdena'] = np.zeros(tdata['edena'].shape)
        tdata['tdenm'] = np.zeros(tdata['edenm'].shape)
        tdata['tdeni'] = np.zeros(tdata['edeni'].shape)

        # tden*: ratio computed only where pden* > 0, zero elsewhere
        #        division by 1.5 according to t*b2 computed in eirmod_extrab25.F90

        here = tdata['pdena'] > 0
        tdata['tdena'][here]  = tdata['edena'][here] / tdata['pdena'][here] / eV / 1.5 # eV
        here = tdata['pdenm'] > 0
        tdata['tdenm'][here]  = tdata['edenm'][here] / tdata['pdenm'][here] / eV / 1.5 # eV
        here = tdata['pdeni'] > 0
        tdata['tdeni'][here]  = tdata['edeni'][here] / tdata['pdeni'][here] / eV / 1.5 # eV

        ########################

        # total neutral quantities
        #
        # CAUTION.
        #
        # - Di-atomic molecules only
        # - validity of t = e / n for multiple species?

        try:
            tdata['pdenn'] = tdata['pdena'].sum(axis = 2)
            tdata['edenn'] = tdata['edena'].sum(axis = 2)
        except:
            tdata['pdenn'] = tdata['pdena'].copy()
            tdata['edenn'] = tdata['edena'].copy()

        try:
            tdata['pdenn'] += 2 * tdata['pdenm'].sum(axis = 2)
            tdata['edenn'] += 2 * tdata['edenm'].sum(axis = 2)
        except:
            tdata['pdenn'] += 2 * tdata['pdenm'].copy()
            tdata['edenn'] += 2 * tdata['edenm'].copy()

        here = tdata['pdenn'] > 0
        tdata['tdenn'] = tdata['edenn'][here] / tdata['pdenn'][here] / eV / 1.5

        ########################

        tdata['vxdena'] = read_ft44_rfield(fid,ver,'vxdena',[ntri,natm])*1e1 # kg s^{-1} m^{-2}
        tdata['vxdenm'] = read_ft44_rfield(fid,ver,'vxdenm',[ntri,nmol])*1e1
        tdata['vxdeni'] = read_ft44_rfield(fid,ver,'vxdeni',[ntri,nion])*1e1

        tdata['vydena'] = read_ft44_rfield(fid,ver,'vydena',[ntri,natm])*1e1 # kg s^{-1} m^{-2}
        tdata['vydenm'] = read_ft44_rfield(fid,ver,'vydenm',[ntri,nmol])*1e1
        tdata['vydeni'] = read_ft44_rfield(fid,ver,'vydeni',[ntri,nion])*1e1

        tdata['vzdena'] = read_ft44_rfield(fid,ver,'vzdena',[ntri,natm])*1e1 # kg s^{-1} m^{-2}
        tdata['vzdenm'] = read_ft44_rfield(fid,ver,'vzdenm',[ntri,nmol])*1e1
        tdata['vzdeni'] = read_ft44_rfield(fid,ver,'vzdeni',[ntri,nion])*1e1


        ## Close file

        fid.close()

        if save is True:
            pickle.dump(clean(tdata), open(os.path.join(where, "fort.46.pkl"), "wb"))

    return clean(tdata)
