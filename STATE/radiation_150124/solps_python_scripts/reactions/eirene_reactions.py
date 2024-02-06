import os
import numpy as np
from solps_python_scripts.reactions.read_reaction import read_reaction

def eirene_reactions(where = "."):

	file = open(os.path.join(where, 'input.dat'), 'r')

	coef = {}

	while True:

		line = file.readline()

		if '*** 4. Data for species and atomic physics module' in line:

			line = file.readline()
			num_reactions = int(file.readline())

			for i in range(num_reactions):

				line = file.readline().split()

				num = line[0]
				coef[num] = {}
				coef[num]['database'] = line[1].swapcase()
				coef[num]['group'] = line[2]
				coef[num]['reaction'] = line[3]
				coef[num]['type'] = line[4]
				try:
					coef[num]['data'] = read_reaction(database = line[1].swapcase(),
													  group    = line[2],
													  reaction = line[3])
				except:
					coef[num]['data'] = {}
					print('  NOT found:' + ' '.join(line))

			break
	
	return coef