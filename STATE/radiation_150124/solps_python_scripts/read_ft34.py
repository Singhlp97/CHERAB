import os
import numpy as np

def read_ft34(where = None):
   
   # cells = read_ft34(file)
   #
   # Read fort.34-files (nodes composing each triangle). 
   #
   #

   # Author: Wouter Dekeyser
   # November 2016
   #
   # Re-writte in python by: Matteo Moscheni
   # E-mail: matteo.moscheni@tokamakenergy.co.uk
   # February 2022

   if os.path.exists(os.path.join(where, "fort.34")) is True:
      fid = open(os.path.join(where, "fort.34"), "r")
   elif os.path.exists(os.path.join(where, "../baserun/fort.34")) is True:
      fid = open(os.path.join(where, "../baserun/fort.34"), "r")
   else: raise ValueError('fort.34 can NOT be found :(')

   ## Read data

   # number of triangels
   ntria = int(fid.readline().split()[0])

   cells = np.zeros((ntria,3), dtype = np.int32)

   for i in range(ntria):
      line = fid.readline().split()
      for j in range(1, len(line)):
         cells[i, j - 1] = int(int(line[j]) - 1)

   # close file
   fid.close()

   return cells
