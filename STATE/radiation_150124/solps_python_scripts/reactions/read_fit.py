import numpy as np

####################################################################

def read_1pFit(index = None, lines = None):

	# Read from AMJUEL database a reaction of the H.1 group

	target = r'\begin{small}\begin{verbatim}'

	line = lines[index]

	while target not in line and index < len(lines)-1:
		index += 1
		line = lines[index]

	if index == len(lines) - 1:
		raise ValueError('Not found...')

	# skip empty or flag line, if any
	
	if lines[index+1] == '' or lines[index+1].split()[0] == 'fit-flag':
		index += 2
	else: index += 1

	# try/except to accomodate AMMONX with 2 coefficients only

	try:

		# read following 3 lines of 3 coefficients

		coef = np.zeros((9))
		for i in range(3):
			line = lines[index + i]
			line = line.replace('D', 'E').split()
			for j in range(3):
				coef[i * 3 + j] = float(line[j * 2 + 1])

	except:

		# read following 1 line of 2 coefficients (AMMONX)
		# with the other left to 0.0

		coef = np.zeros((9))
		line = lines[index]
		line = line.replace('D', 'E').split()
		for j in range(2):
			coef[j] = float(line[j * 2 + 1])

	return coef

####################################################################             

def read_2pFit(index = None, lines = None):

	# Read from AMJUEL database a reaction of the H.1 group

	#target = r'\begin{small}\begin{verbatim}'
	target_1 = r'T-Index'
	target_2 = r'T Index'

	line = lines[index]

	while target_1 not in line and target_2 not in line and index < len(lines)-1:
		index += 1
		line = lines[index]

	if index == len(lines) - 1:
		raise ValueError('Not found...')

	# read following lines of coefficients (skip headers)

	# - axis = 0: T
	# - axis = 1: E

	index += 1

	coef = np.zeros((9,9))
	for block in range(3):
		for i in range(9):
			line = lines[index]
			line = line.replace('D', 'E').split()
			for j in range(3):
				coef[i,j + 3 * block] = float(line[j + 1])
			index += 1
		index += 3

	return coef               