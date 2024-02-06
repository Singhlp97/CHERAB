import numpy as np
import scipy.special as sp

#############################

# constants: I.S. units

e = 1.602E-19
eps0 = 8.854E-12
me = 9.109E-31 
c = 3E+08
h = 6.626E-34
hbar = h / 2 / np.pi

#############################

# plasma quantities

ne = 1E+19 # ????
ni = np.array([1E+19])

if len(ni) == 1 and ne != ni[0]:
    raise ValueError('ni = ne required for hydrogenic plasmas!')

te = 2E+03 * e # ????
Zi = np.array([1])
Z2ni = (Zi**2 * ni).sum()

def wp(ne):
    return np.sqrt(ne * e**2 / eps0 / me)

#############################

# camera quantities

# Ha = 656 nm
lHa =     656E-09
lHa_min = 650E-09
lHa_max = 670E-09

def to_frequency(l):
    return 2 * np.pi * c / l 

wHa = to_frequency(lHa)
wHa_max = to_frequency(lHa_min)
wHa_min = to_frequency(lHa_max)

#############################

# mathematics [https://en.wikipedia.org/wiki/Bremsstrahlung]

def C1():
    return 8 / 3 * np.sqrt(2 / np.pi) * (e**2 / 4 / np.pi / eps0)**3 / (me * c**2)**(1.5)
def C2(w, ne):
    return np.sqrt(1 - (wp(ne) / w)**2)
def C3(ne, te, Z2ni):
    return Z2ni * ne / np.sqrt(te) 

def y(w):
    return 0.5 * (hbar * w / te)**2

# E1 asympthotically approximated when y << 1

def E1(y):
    return sp.expi(y)
    # return - np.log(y * np.exp(0.577)) # valid only if y << 1

#############################

# computation

def dPdw(w):
    return C1() * C2(w, ne) * C3(ne, te, Z2ni) * E1(y(w))

w = np.linspace(wHa_min, wHa_max, int(1E+02))
dw = np.abs(w[:-1] - w[1:])

PbremHa = (dPdw(w[:-1]) * dw).sum() # per unit volume

# - minimum w at wp(ne): if below, radiation can NOT make it through
# - maximum at 2E+18, where y ~ 0.1 (i.e. y << 1 approximation limit)

w = np.logspace(np.log10(wp(ne)), np.log10(2E+18), int(1E+04))
dw = np.abs(w[:-1] - w[1:])

PbremTot = (dPdw(w[:-1]) * dw).sum()

#############################

print()
print('Bremsstrahlung fraction in Ha spectrum: {:.2G} %'.format(PbremHa / PbremTot * 1E+02))
print()



