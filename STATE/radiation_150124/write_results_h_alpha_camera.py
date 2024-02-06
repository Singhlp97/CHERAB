parula = [[0.2081, 0.1663, 0.5292], 
                [0.2116238095, 0.1897809524, 0.5776761905], 
                [0.212252381, 0.2137714286, 0.6269714286], 
                [0.2081, 0.2386, 0.6770857143], 
                [0.1959047619, 0.2644571429, 0.7279], 
                [0.1707285714, 0.2919380952, 0.779247619], 
                [0.1252714286, 0.3242428571, 0.8302714286], 
                [0.0591333333, 0.3598333333, 0.8683333333], 
                [0.0116952381, 0.3875095238, 0.8819571429], 
                [0.0059571429, 0.4086142857, 0.8828428571], 
                [0.0165142857, 0.4266, 0.8786333333], 
                [0.032852381, 0.4430428571, 0.8719571429], 
                [0.0498142857, 0.4585714286, 0.8640571429], 
                [0.0629333333, 0.4736904762, 0.8554380952], 
                [0.0722666667, 0.4886666667, 0.8467], 
                [0.0779428571, 0.5039857143, 0.8383714286], 
                [0.079347619, 0.5200238095, 0.8311809524], 
                [0.0749428571, 0.5375428571, 0.8262714286], 
                [0.0640571429, 0.5569857143, 0.8239571429], 
                [0.0487714286, 0.5772238095, 0.8228285714], 
                [0.0343428571, 0.5965809524, 0.819852381], 
                [0.0265, 0.6137, 0.8135], 
                [0.0238904762, 0.6286619048, 0.8037619048], 
                [0.0230904762, 0.6417857143, 0.7912666667], 
                [0.0227714286, 0.6534857143, 0.7767571429], 
                [0.0266619048, 0.6641952381, 0.7607190476], 
                [0.0383714286, 0.6742714286, 0.743552381], 
                [0.0589714286, 0.6837571429, 0.7253857143], 
                [0.0843, 0.6928333333, 0.7061666667], 
                [0.1132952381, 0.7015, 0.6858571429], 
                [0.1452714286, 0.7097571429, 0.6646285714], 
                [0.1801333333, 0.7176571429, 0.6424333333], 
                [0.2178285714, 0.7250428571, 0.6192619048], 
                [0.2586428571, 0.7317142857, 0.5954285714], 
                [0.3021714286, 0.7376047619, 0.5711857143], 
                [0.3481666667, 0.7424333333, 0.5472666667], 
                [0.3952571429, 0.7459, 0.5244428571], 
                [0.4420095238, 0.7480809524, 0.5033142857], 
                [0.4871238095, 0.7490619048, 0.4839761905], 
                [0.5300285714, 0.7491142857, 0.4661142857], 
                [0.5708571429, 0.7485190476, 0.4493904762],
                [0.609852381, 0.7473142857, 0.4336857143], 
                [0.6473, 0.7456, 0.4188], 
                [0.6834190476, 0.7434761905, 0.4044333333], 
                [0.7184095238, 0.7411333333, 0.3904761905], 
                [0.7524857143, 0.7384, 0.3768142857], 
                [0.7858428571, 0.7355666667, 0.3632714286], 
                [0.8185047619, 0.7327333333, 0.3497904762], 
                [0.8506571429, 0.7299, 0.3360285714], 
                [0.8824333333, 0.7274333333, 0.3217], 
                [0.9139333333, 0.7257857143, 0.3062761905], 
                [0.9449571429, 0.7261142857, 0.2886428571], 
                [0.9738952381, 0.7313952381, 0.266647619], 
                [0.9937714286, 0.7454571429, 0.240347619], 
                [0.9990428571, 0.7653142857, 0.2164142857], 
                [0.9955333333, 0.7860571429, 0.196652381], 
                [0.988, 0.8066, 0.1793666667], 
                [0.9788571429, 0.8271428571, 0.1633142857], 
                [0.9697, 0.8481380952, 0.147452381], 
                [0.9625857143, 0.8705142857, 0.1309], 
                [0.9588714286, 0.8949, 0.1132428571], 
                [0.9598238095, 0.9218333333, 0.0948380952], 
                [0.9661, 0.9514428571, 0.0755333333], 
                [0.9763, 0.9831, 0.0538]]

import os
import csv
import math
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap

from raysect.primitive.mesh.vtk import export_vtk

####################################################################################################

def plot(cfg = None, name = None, x = None, y = None, data = None, title = None):

    basename_dir = os.path.join(cfg['output_directory_extended'], '..')

    # colormap: greys reversed (from black to white)
    cmap = plt.cm.get_cmap('Greys')
    reversed_cmap = cmap.reversed()

    fig, ax = plt.subplots()
    c = ax.pcolormesh(x, y[::-1], data, cmap = reversed_cmap, vmin = data.min(), vmax = data.max())
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    #plt.axis('equal')
    plt.axis('off')
    fig.colorbar(c, ax = ax)  
    ax.set_title(title)
    plt.savefig(os.path.join(basename_dir, name + '.png'))

    return

####################################################################################################

def write_figures(cfg = None, name = None, x = None, y = None, data = None, data_err = None):

    # heat load lin

    plot(cfg = cfg, name = name + "_lin", x = x, y = y, data = data, title = 'H_alpha camera view: lin scale')

    # heat load log

    data_log = np.log10(data + data[data > 0].min() * 1e-1)
    plot(cfg = cfg, name = name + "_log", x = x, y = y, data = data_log, title = 'H_alpha camera view: log scale')
    
    # heat load error lin

    data[data == 0] = np.nan
    err = np.abs(data_err / data)
    err[np.isnan(err)] = 0

    plot(cfg = cfg, name = name + "_error_lin", x = x, y = y, data = err, title = 'H_alpha camera view error: lin scale')

    # heat load error log

    err = np.log10(err + err[err > 0].min() / 10)
    plot(cfg = cfg, name = name + "_error_log", x = x, y = y, data = err, title = 'H_alpha camera view error: log scale')

####################################################################################################
####################################################################################################
####################################################################################################

def write_results(cfg = None, name = None, camera = None, pipeline = None, pixel_origins = None):

    side = 2 * np.tan(np.pi / 180 * 0.5 * \
           cfg['raytracing']['observer']['specs']['fov']) * \
           cfg['raytracing']['observer']['surface_offset']

    if pixel_origins is None:

        nx = cfg['raytracing']['observer']['specs']['nx_pixels']
        ny = cfg['raytracing']['observer']['specs']['ny_pixels']

        x, y = np.meshgrid(np.linspace(0, side, nx), np.linspace(0, side, ny))
        (y, x) = (x.T, y.T)

    else:

        x = [0]
        xcum = 0
        y = [0]
        ycum = 0

        for i in range(np.size(pixel_origins,0) - 1):
            
            p1 = pixel_origins[i][0]
            p2 = pixel_origins[i+1][0]
            lx12 = p1.distance_to(p2)
            x.append(lx12 + xcum)
            xcum += lx12

        for i in range(np.size(pixel_origins,1) - 1):

            p1 = pixel_origins[0][i]
            p2 = pixel_origins[0][i+1]
            ly12 = p1.distance_to(p2)
            y.append(ly12 + ycum)
            ycum += ly12

        print(len(x))

        ############

        (x,y) = (np.array(x), np.array(y))

        plt.figure()
        for i in range(len(x)):
            plt.plot([y[0], y[-1]], [x[i], x[i]], 'm-')
        for i in range(len(y)):
            plt.plot([y[i], y[i]], [x[0], x[-1]], 'm-')
        plt.title('pixel grid')
        plt.axis('off')
        plt.savefig(os.path.join(cfg['output_directory_extended'], 'pixel_grid.png'))
        # plt.show()

        ############

        x, y = np.meshgrid(x, y)
        (x, y) = (x.T, y.T)

    frame = pipeline.frame

    # from W/m2/str to W/m2
    heatLoad = frame.mean
    heatLoad_errors = frame.errors()

    heatLoad[np.isnan(heatLoad)] = 0
    heatLoad_errors[np.isnan(heatLoad_errors)] = 0
    # CAUTION: convention of the pixel ordering is NOT physical-space like
    heatLoad = heatLoad[:, ::-1]
    heatLoad_errors = heatLoad_errors[:, ::-1]

    basename_dir = os.path.join(cfg['output_directory_extended'], '..')

    # pixel_grid
    if pixel_origins is not None:
        pickle.dump((x,y), open(os.path.join(basename_dir, 'input', 'pixel_grid.pkl'), 'wb'))

    # heat load view
    with open(os.path.join(basename_dir, name + '.csv'), 'w+', newline = '') as f:
        writer = csv.writer(f)
        for i in range(frame.nx):
            writer.writerow((heatLoad[i,:]))
    with open(os.path.join(basename_dir, name + '_errors.csv'), 'w+', newline = '') as f:
        writer = csv.writer(f)
        for i in range(frame.nx):
            writer.writerow((heatLoad_errors[i,:]))

    # total power received
    #with open(os.path.join(basename_dir, 'Ptot_' + name + '.csv'), 'w+', newline = '') as f:
    #    writer = csv.writer(f)
    #    writer.writerow([np.sum(heatLoad * side / nx * side / ny)])

    #write_figures(cfg = cfg, name = name, data = heatLoad, data_err = heatLoad_errors,
    #              x = np.linspace(0, nx, nx+1), y = np.linspace(0, ny, ny+1))
    write_figures(cfg = cfg, name = name, data = heatLoad, data_err = heatLoad_errors,
                  x = x, y = y)

    lux_file = 'illuminance'
    radius_oblo = 0.5
    side = 2 * np.tan(np.pi / 180 * 0.5 * cfg["raytracing"]["observer"]["specs"]["fov"]) * 0.9829523234899655

    x_shift = x - (x[0,-1]/2 + 0.03528361056281585)
    y_shift = y - (y[-1,-1]/2 + 0.08541462678416768)

    # for i in range(len(y)):
    #     for j in range(len(x)):
    #         if np.sqrt(x_shift[i][j]**2+y_shift[i][j]**2) > radius_oblo:
    #             heatLoad[i][j] = 0 

    parula_map = LinearSegmentedColormap.from_list('parula', parula)

    conversion_ptoW= 3e-19
    W_m2 = heatLoad*conversion_ptoW
    conversionWtolumen = 1.019*math.exp(-285.4*pow((0.656-0.559),2))
    lumen_m2 = W_m2*conversionWtolumen*683*2*math.pi # for 656 nm 

    lux = lumen_m2
    lux[np.isnan(lux)] = 0
    
    with open(os.path.join(basename_dir, lux_file + ".csv"), "w+", newline = "") as f:
        writer = csv.writer(f)
        for i in range(1520):
            writer.writerow((lux[i,:]))

    lux[lux<1] = 1
    # plt.figure(199)
    # plt.imshow(lux, norm=LogNorm(vmin=1, vmax=lux.max()), cmap=parula_map)

    fig, ax = plt.subplots()
    c = ax.pcolormesh(x, y[::-1], lux, cmap = parula_map, norm=LogNorm(vmin=1, vmax=lux.max()))
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    #plt.axis('equal')
    plt.axis('off')
    fig.colorbar(c, ax = ax)  
    ax.set_title('Illuminance [lux]')
    plt.savefig(os.path.join(basename_dir, lux_file + '.png'))

    print('Max illuminance = ', lux.max())



    # plot(cfg = cfg, name = name + "_error_log", x = x, y = y, data = err, title = 'H_alpha camera view error: log scale')
    # cb = plt.colorbar()
    # cb.set_label(label='Illuminance [lux]', size=20)
    # plt.axis('off')
    # cb.ax.tick_params(labelsize=20)
    # #plt.title('Illuminance [lux]', fontsize=15)
    # plt.savefig(os.path.join(basename_dir, lux_file + ".png"))