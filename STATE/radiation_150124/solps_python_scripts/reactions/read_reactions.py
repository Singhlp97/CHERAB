import os
import pickle
import warnings
import numpy as np

from solps_python_scripts.reactions.read_fit import read_1pFit, read_2pFit

#########################################################################################

def read_one_reaction(database = None, group = None, reaction = None, dictionary = None):

    TREETOP = '/home/lsingh/runs_20220929_thinplasma/solps_python_scripts/reactions'

    fid = open(os.path.join(TREETOP, database), 'r') 
    lines = [line for line in fid]

    group_label = ''.join(['\\section{', group])
    reaction_label = ''.join(['Reaction ', reaction])

    group_found = False
    coef = []

    for index, line in enumerate(lines):

        if group_label in line: group_found = True

        if group_found is True:

            # end of the right group, search aborted, reaction not found

            # search for reaction

            if reaction_label in line:

                # any suffix (e.g. 2.1.5FJ) is correctly included/excluded

                if line.split()[1] == reaction:

                    # try/except to avoid error for AMMONX not having $-split formulae
                    try:    equation = line.split('$')[1].replace('\\ge', '>=')
                    except: equation = '-'

                    if dictionary is not None:
                        dictionary['equation']  = ''.join(['$', equation, '$'])
                        for isp in range(2):
                            # different cases because of AMMONX convention
                            if 'R-H' in reaction: species = reaction.split('-')[isp+1]
                            else:                 species = equation.split('\\rightarrow')[0].split(' + ')[isp].strip()
                            # change name according to fort.44 notation, i.e. remove any "_" or "^"
                            # for symbol in ['_', '^']: species = ''.join(species.split(symbol))
                            dictionary['species_' + str(isp)] = species

                    if   group == 'H.0':
                        coef = read_1pFit(index = index, lines = lines)
                    elif group == 'H.1':
                        coef = read_1pFit(index = index, lines = lines)
                    elif group == 'H.2' and database != 'AMMONX':
                        coef = read_1pFit(index = index, lines = lines)
                    elif group == 'H.3':
                        coef = read_2pFit(index = index, lines = lines)
                    elif group == 'H.4':
                        coef = read_2pFit(index = index, lines = lines)
                    elif group == 'H.8':
                        coef = read_1pFit(index = index, lines = lines)
                    elif group == 'H.9':
                        coef = read_2pFit(index = index, lines = lines)
                    elif group == 'H.10':
                        coef = read_2pFit(index = index, lines = lines)
                    elif group == 'H.11':
                        coef = read_1pFit(index = index, lines = lines)
                    elif group == 'H.12':
                        coef = read_2pFit(index = index, lines = lines)
                    elif database == 'AMMONX':
                        coef = read_1pFit(index = index, lines = lines)

                    # read range of validity (if any)

                    # CAUTION.
                    #
                    # - 4 rows at most after coefficient block
                    # - '=' and ('MIN' or 'MAX') must be present
                    # - density in [cm^{-3}] => [m^{-3}]

                    break

            if ''.join(['\\section{', group.split('.')[0] + '.' + str(int(group.split('.')[1]) + 1)]) in line:
                warnings.warn('Not found: {} {} {}'.format(database, group, reaction))
                break

    if dictionary is not None: dictionary['fit_coefficients'] = coef

    return dictionary or coef
                

#########################################################################################

def read_all_reactions(where = ".", verbose = True, save = True):

    try:

        coef = pickle.load(open(os.path.join(where, 'reactions.pkl'), 'rb'))
        if verbose is True:
            print('read_all_reactions - reactions.pkl')
            for ID in coef.keys():
                try:
                    if coef[ID]['fit_coefficients'] == {}:
                        label = ' '.join([coef[ID]['number'],
                                          coef[ID]['database'],
                                          coef[ID]['group'],
                                          coef[ID]['reaction'],
                                          coef[ID]['type']])
                        print('  NOT found: ' + label)
                except: pass
        return coef 

    except:

        if verbose is True: print('read_all_reactions - input.dat')

        # file = open(os.path.join(where, 'input.dat'), 'r')
        # lines = file.read()
        # num_lines = len(lines)
        # file.close()
        file = open(os.path.join(where, 'input.dat'), 'r')

        i_line = 0
        coef = {}

        while True:

            line = file.readline()
            i_line += 1

            if '*** 4. Data for species and atomic physics module' in line:

                line = file.readline()
                num_reactions = int(file.readline())
                num = 0
                ID  = str(0)

                while int(num) < num_reactions:

                    # try/except to avoid dummy line (e.g. N    1)

                    try:

                        line = file.readline().split()

                        num_old = num
                        num = line[0]

                        coef[ID]              = {}
                        coef[ID]['number']    = line[0]
                        coef[ID]['database']  = line[1]

                        # try/except to avoid error when H.10 x.y.z
                        # spelled withotuh blank spaces (H.10x.y.z)
                        try:
                            coef[ID]['group']     = line[2]
                            coef[ID]['reaction']  = line[3]
                            coef[ID]['type']      = line[4]
                            coef[ID]['ispecies_0'] = int(line[5])
                            coef[ID]['ispecies_1'] = int(line[6])
                        except:
                            coef[ID]['group']     = line[2][:4]
                            coef[ID]['reaction']  = line[2][4:]
                            coef[ID]['type']      = line[3]
                            coef[ID]['ispecies_0'] = int(line[4])
                            coef[ID]['ispecies_1'] = int(line[5])

                        # try/except to avoid AMMONX/ADAS databases

                        try:
                            coef[ID] = read_one_reaction(database   = line[1],
                                                         group      = line[2],
                                                         reaction   = line[3],
                                                         dictionary = coef[ID])
                        except:
                            coef[ID]['fit_coefficients'] = {}
                            print('  NOT found: ' + ' '.join(line[:5]))

                        ID = str(int(ID) + 1)

                    except: num = num_old

                break

        if save is True:
            pickle.dump(coef, open(os.path.join(where, 'reactions.pkl'), 'wb'))
        
    return coef
