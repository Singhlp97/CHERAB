import os
import numpy as np

def read_ft35(where = None):

   # necke = read_ft35(file)
   #
   # Read fort.35-files (triangle data). 
   #
   #

   # Author: Wouter Dekeyser
   # November 2016
   #
   # Re-writte in python by: Matteo Moscheni
   # E-mail: matteo.moscheni@tokamakenergy.co.uk
   # February 2022

   if os.path.exists(os.path.join(where, "fort.35")) is True:
      fid = open(os.path.join(where, "fort.35"), "r")
   elif os.path.exists(os.path.join(where, "../baserun/fort.35")) is True:
      fid = open(os.path.join(where, "../baserun/fort.35"), "r")
   else: raise ValueError('fort.35 can NOT be found :(')

   ## Read data

   # number of triangels
   ntria = int(fid.readline().split()[0])

   links = {}

   nghbr = np.zeros((ntria,3))
   side  = np.zeros((ntria,3))
   cont  = np.zeros((ntria,3))
   ixiy  = np.zeros((ntria,2))

   for i in range(ntria):
      data = np.array(fid.readline().split())
      nghbr[i,:] = [int(data[1]), int(data[4]), int(data[7])]
      side[i,:]  = [int(data[2]), int(data[5]), int(data[8])]
      cont[i,:]  = [int(data[3]), int(data[6]), int(data[9])]
      ixiy[i,:]  = [int(data[10]), int(data[11])]

   links['nghbr'] = nghbr
   links['side']  = side
   links['cont']  = cont
   links['ixiy']  = ixiy

   # close file
   fid.close()

   return links