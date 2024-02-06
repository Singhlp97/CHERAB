import numpy as np
import matplotlib.pyplot as plt

######################################
# absorption cross section: n = 2 -> 3

# bibliography:
# [1] "Accurate atomic transition probabilities for hydrogen, helium and lithium",
#      W.L. Wiese, J.R. Fuhr, 2009.
# [2] "Line strengths, A-factors and absorption cross sections for...",
#      O. Axner, et al., 2004.
# [3] "Inelastic collision rates of trapped metastable hydrogen",
#      D. Landhuis, et al., 2008.
# [4] "On inelastic hydrogen atom collisions in stellar atmosphere",
#      P.S. Barklem, et al., 2011.

# DATA     from ref. [1]
# THEORY   from ref. [2]
# NOTATION from ref. [2]
#
# Ha = absorption from n=2 to n=3

# SI units

e    = 1.60217E-19
eps0 = 8.85419E-12
me   = 9.10938E-31
c    = 2.99792E+08
h    = 6.62607E-34

#l_Ha = 656E-09   # [m]
#nu_Ha = c / l_Ha # [Hz]

# oscillator strength
f_Ha  = 6.41080E-01
s0_tilde = e**2 / 4 / eps0 / me / c

# from [1]
A_21 = 4.6986E+00 * 1E+08 # = A_12
A_31 = 5.5751E-01 * 1E+08 # = A_13
A_32 = 4.4101E-01 * 1E+08 # = A_23

############# INELASTIC PROCESSES ##############

n_H0 = 1E+19 # [1/m^3]: very high => conservativeness
T_eV = np.logspace(0,2,250)

# from [3]: upper bound => conservativeness ([m^3/s])
# transition: n=2 => n=1
sv_21_inel = (1.8 + 1.8) * 1E-15

# from [4]: Darwin function => overestimate
#           => conservativeness ([m^3/s])
# transition: n=3 => n={2,1}
#
# Hydrogen-to-Hydrogen inelastic scattering
# => mA = mH => simplifications

# E_ul = 2.0        # Na 3p-3s
# f_lu = 0.1051E+01 # Na 3p-3s

def sv_Darwin_inel(n_u = 0, n_l = 0, f_ul = 0, T_eV = None):

	a0  = 5.29177E-11 # [m] - Bohr radius
	E1H = 13.6        # [eV] - H ionisation potential
	me  = 9.10938E-31 # [kg] - electron mass (rest)
	mp  = 1.67262E-27 # [kg] - proton mass (rest)
	mH  = mp + me     # [kg] - hydrogen mass 
	e   = 1.60217E-19 # [C]  - electron charge
	mu  = mp * me / (mp + me)

	T_J = T_eV * e # [J] = kB * T_K

	def E(n):
		return -E1H / n**2 # [eV]

	E_ul = np.abs(E(n_u) - E(n_l)) # [eV]
	w_ul = E_ul / T_eV             # [eV]

	def psi(w_ul):
		return np.exp(-w_ul) * (1 + 2 / w_ul)

	return 16 * np.pi * a0**2 * (E1H / E_ul) * f_ul * \
		   me * (2 * mH) / mH / (mH + me) *		  \
		   np.sqrt(2 * T_J / np.pi / mu) * psi(w_ul)

# 3 -> 2 transition
f_32 = 6.4108E-01
sv_32_inel = sv_Darwin_inel(n_u = 3, n_l = 2, f_ul = f_32, T_eV = T_eV)

# 3 -> 1 transition
f_31 = 7.9142E-02
sv_31_inel = sv_Darwin_inel(n_u = 3, n_l = 1, f_ul = f_31, T_eV = T_eV)

plt.figure()
plt.loglog(T_eV, n_H0 * sv_32_inel, label = "3 -> 2")
plt.loglog(T_eV, n_H0 * sv_31_inel, label = "3 -> 1")
plt.xlabel('$T_{H_0}\;[eV]$')
plt.ylabel(r'$n\cdot \langle\sigma v\rangle\;[s^{-1}]$')
plt.ylim(top = 1E+07)
plt.title('Inelastic rates from $n=3$')
plt.legend(loc = 'lower right')
plt.grid(True)

############# ELASTIC PROCESSES ##############

print()
print("WARNING. gamma @resonance: neglecting elastic processes => NON-conservative")
print()

T_ref = 1E+1 # [eV]
sv_32_inel = sv_Darwin_inel(n_u = 3, n_l = 2, f_ul = f_32, T_eV = T_ref)
sv_31_inel = sv_Darwin_inel(n_u = 3, n_l = 1, f_ul = f_31, T_eV = T_ref)

gamma_at_resonance = A_21 + (A_31 + A_32) + n_H0 * (sv_21_inel + sv_31_inel + sv_32_inel)

s_absorption_Ha = (4 / gamma_at_resonance) * s0_tilde * f_Ha

print("Halpha absorption cross section: {:.4G} m^2 @T_H0={:.3G}eV".format(s_absorption_Ha, T_ref))
print()
print("NEGLECTING ELASTIC SCATTERING: 1E4 vs 1E8 s^{-1}")
print()

plt.show()