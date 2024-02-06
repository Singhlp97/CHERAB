# homemade linear fit of fig. 1 and 3 in:
# https://www.researchgate.net/publication/240388531_Temperature_dependence_of_the_emissivity_of_transition_metals
#
# total hemispherical emissivities considered

def coefficients(t1, t2, e1, e2):
	m = (e2 - e1) / (t2 - t1)
	q = e1 - m * t1
	return m, q

def EmissivityMolybdenum(T):

	# temperature = 500 K
	# emissivity = 0.05
	t1 = 5.00E+02
	e1 = 5.00E-02

	# temperature = 2300 K
	# emissivity = 0.25
	t2 = 2.30E+03
	e2 = 2.50E-01

	# straight line coefficients
	m, q = coefficients(t1, t2, e1, e2)

	return m * T + q

def EmissivityTungsten(T):

	# temperature = 300 K
	# emissivity = 0.05
	t1 = 3.00E+02
	e1 = 5.00E-02

	# temperature = 2300 K
	# emissivity = 0.25
	t2 = 2.75E+03
	e2 = 3.00E-01

	# straight line coefficients
	m, q = coefficients(t1, t2, e1, e2)

	return m * T + q