import os
import numpy as np
import math

#######################################################################################################################

#location = 'radiation_load/output/Bypass_leakage/total_radiation/2.step=True.abs=False.astra=False.b2=True.eirene=False.AxisymmetricMapper/'
# location = 'radiation_load/output/Dpuff_Npuff_PFR.run.1.BGK=OFF.BCCON=8.2.5E20.Dpuff=1E21.Npuff=6E20.ANcoeff.LJ-VA/total_radiation/1.step=True.abs=False.astra=False.b2=True.eirene=False.extra=False.AxisymmetricMapper/'
# location = 'radiation_load/output/RF.mirror.VS_N_puff=6.00E19_Pe=2.27MW_Pi=1.49MW_sequence.run.13/total_radiation/1.step=True.abs=False.astra=False.b2=True.eirene=False.extra=False.AxisymmetricMapper/'
location = 'radiation_load/output/Fabio/total_radiation/1.step=True.abs=False.scat=False.astra=False.b2=True.eirene=False.extra=False.AxisymmetricMapper'

run_files = os.listdir(location)
rel_errors = []
for file in run_files:
    filename, ext = os.path.splitext(file)
    if ext == '.csv' and filename == 'camera_shift':
        print("importing", file)
        powers = np.loadtxt(open(os.path.join(location + '/', filename + '.csv'), "r"), delimiter=",", skiprows=1)
        try:
            errors = np.loadtxt(open(os.path.join(location + '/', filename + '_errors.csv'), "r"), delimiter=",", skiprows=1)
            power_m2_errors = errors[:]
            powers_m2 = powers[:]
            num_detectors = np.size(powers[:,-1])
            for i in range(num_detectors):
                for j in range(num_detectors):
                    if powers_m2[i][j] > 0:
                        rel_errors += [power_m2_errors[i][j] / powers_m2[i][j]]
        except:
            powers_m2 = powers[:,1]
            power_m2_errors = powers[:,2]
            num_detectors = np.size(powers[:,-1])
            for i in range(num_detectors):
                if powers_m2[i] > 0:
                    rel_errors += [power_m2_errors[i] / powers_m2[i]]
        
AVG_error = np.mean(rel_errors)
print('\n\nAverage error: {:.8G} %'.format(AVG_error * 100))
print()
