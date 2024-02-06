import os

def read_b2mn(where = "."):

	b2mn = {}

	f = open(os.path.join(where, 'b2mn.dat'))

	for line in f:
		line = line.split()
		try:    b2mn[line[0].split("'")[1]] = float(line[1].split("'")[1])
		except: pass

	if 'b2mwti_jxi' not in b2mn.keys(): b2mn['b2mwti_jxi'] = 0

	return b2mn