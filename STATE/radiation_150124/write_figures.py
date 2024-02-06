# %% 
import os
import json
import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
from raysect.primitive.mesh.vtk import export_vtk

parser = argparse.ArgumentParser('Calculate power deposition.')
parser.add_argument('configFile', type = str, help = "A JSON format task configuration file.")
args = parser.parse_args()
with open(args.configFile, 'r') as f:
    cfg = json.load(f)

basename_dir = './output/H_alpha/hal/'
basename_file = 'camera_shift'

heatLoad = np.loadtxt(basename_dir + basename_file + '.csv', delimiter=',')
heatLoad_errors = np.loadtxt(basename_dir + basename_file + '_errors.csv', delimiter=',')

side = 2 * np.tan(np.pi / 180 * 0.5 * cfg["camera"]["observation_angle"]) * 0.9829523234899655

nx = cfg["camera"]["nx_pixels"]
ny = cfg["camera"]["ny_pixels"]
x, y = np.meshgrid(np.linspace(0, side, nx), np.linspace(0, side, ny))

# colormap: greys reversed (from black to white)
cmap = plt.cm.get_cmap('Greys')
reversed_cmap = cmap.reversed()

fig, ax = plt.subplots()

c = ax.pcolormesh(x, y[::-1], heatLoad, cmap = reversed_cmap, vmin = heatLoad.min(), vmax = heatLoad.max())
ax.axis([x.min(), x.max(), y.min(), y.max()])
fig.colorbar(c, ax = ax)  
ax.set_title('H-alpa camera view - LIN')
plt.savefig(os.path.join(basename_dir, basename_file + ".png"))
plt.close()

heatLoad_log = np.log10(heatLoad + heatLoad[heatLoad > 0].min() * 1e-1)

fig, ax = plt.subplots()

c = ax.pcolormesh(x, y[::-1], heatLoad_log, cmap = reversed_cmap, vmin = heatLoad_log.min(), vmax = heatLoad_log.max())
ax.axis([x.min(), x.max(), y.min(), y.max()])
fig.colorbar(c, ax = ax)  
ax.set_title('H-alpa camera view - LOG')
plt.savefig(os.path.join(basename_dir, basename_file + "_log.png"))
plt.close()

err = np.abs(heatLoad_errors / heatLoad)
err[np.isnan(err)] = 0

fig, ax = plt.subplots()

c = ax.pcolormesh(x, y[::-1], err, cmap = reversed_cmap, vmin = err.min(), vmax = err.max())
ax.axis([x.min(), x.max(), y.min(), y.max()])
fig.colorbar(c, ax = ax)  
ax.set_title('H-alpa camera view error - LIN')
plt.savefig(os.path.join(basename_dir, basename_file + "_errors_lin.png"))
plt.close()

try:
    err = np.log10(err + err[err > 0].min() / 10)
    fig, ax = plt.subplots()
    c = ax.pcolormesh(x, y[::-1], err, cmap = reversed_cmap, vmin = err.min(), vmax = err.max())
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(c, ax = ax)  
    ax.set_title('H-alpa camera view error - LOG')
    plt.savefig(os.path.join(basename_dir, basename_file + "_errors_log.png"))
    plt.close()

except:
    err = []
# %%

with open(os.path.join(basename_dir, "Ptot_" + basename_file + ".csv"), "w", newline = "") as f:
    writer = csv.writer(f)
    writer.writerow([np.sum(heatLoad * side/nx * side/ny)])
    writer.writerow([np.sum(heatLoad_errors * side/nx * side/ny)])

# %%
plt.show()
# %%

# %%
