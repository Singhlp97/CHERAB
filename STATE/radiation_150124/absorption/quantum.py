import numpy as np
import matplotlib.pyplot as plt


def fun(E, T):
  return np.exp(- np.abs(E + 13.6) / T)  


# def PBsplit(ml, ms):

#   e =  1.6E-19 # [C]
#   hbar = 1.05E-34
#   me = 9E-31

#   gs = 2.0 # [-]
#   Bz = 2.0 # [T]
#   muB= 0.5 * e * hbar / me

#   return np.abs(Bz * muB * (ml + gs * ms) / e)


def E(n):
  return -13.6 / n**2 # [eV]

def g(n):
  return 2 * n**2 # degeneracy

# [-] principal quantum number
n = np.arange(50) + 1

# [eV] H0 temperature
T = np.logspace(0,3,500)

# partition function
Z = np.zeros(T.shape)

for i in range(T.shape[0]):
  for j in range(n.shape[0]):
    Z[i] += fun(E(n[j]),T[i]) * g(n[j])

# dE = np.zeros(n.shape)

# ms = np.array([-0.5, +0.5])

# for i in range(len(n)):
  
#   l = np.arange(n[i])
#   ml = np.arange(-l.max(), l.max()+1)

#   split = []

#   for iml in range(len(ml)):
#     for ims in range(len(ms)):
#       split += [PBsplit(ml[iml], ms[ims])]

#   dE[i] = np.max(split)

plt.figure()

# [-] reference level
refs = np.array([1, 2, 3])
fraction = np.zeros((T.shape[0], refs.shape[0]))

for i in range(len(refs)):

  n = refs[i]
  fraction[:,i] = g(n) * fun(E(n),T) / Z

  plt.loglog(T, fraction[:,i], '-', label = '$H_0(n=' + str(n) + ')$')

plt.xlabel('$T_{H_0}\;[eV]$')
plt.ylabel('$N^{*}/N\;[-]$')
plt.title('Fraction of excited $H_0$')
plt.legend(loc = 'upper right')
plt.grid(True)
plt.savefig("./fraction_excited_atoms.png")

##################################
# n=2 vs. n=3

# bibliography:
# [1] "Accurate atomic transition probabilities for hydrogen, helium and lithium",
#      W.L. Wiese, J.R. Fuhr, 2009.
# [2] "Line strengths, A-factors and absorption cross sections for...",
#      O. Axner, et al., 2004.

A_13 = 5.5751E-01 # [1/s]
A_23 = 4.4101E-01 # [1/s]

fraction_32 = A_23 / (A_13 + A_23)

#g_12 = g(1) + g(2)
#fraction_32 = g(1)/g_12 * A_13 / (g(1)/g_12 * A_13 + g(2)/g_12 * A_23)

N2toN3 = fraction[:,1] / fraction[:,2] # g(2) / g(3) * fun(E(2),T) / fun(E(3),T)
N2toN3 /= fraction_32

plt.figure()
plt.semilogx(T, N2toN3, 'm-')
plt.xlabel('$T_{H_0}\;[eV]$')
plt.ylabel('$N_2\:/\:(A_{32}N_3)\;[-]$')
plt.title(r'$H_{\alpha}$ absorbers ($N_2$) vs. $H_{\alpha}$ emitters ($f_{32}N_3$)')
plt.grid(True)
plt.savefig("./fraction_absorbing_atoms.png")

######################################

plt.show()



