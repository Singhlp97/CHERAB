
# Copyright 2016-2021 Euratom
# Copyright 2016-2021 United Kingdom Atomic Energy Authority
# Copyright 2016-2021 Centro de Investigaciones Energéticas, Medioambientales y Tecnológicas
#
# Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the
# European Commission - subsequent versions of the EUPL (the "Licence");
# You may not use this work except in compliance with the Licence.
# You may obtain a copy of the Licence at:
#
# https://joinup.ec.europa.eu/software/page/eupl5
#
# Unless required by applicable law or agreed to in writing, software distributed
# under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied.
#
# See the Licence for the specific language governing permissions and limitations
# under the Licence.

import numpy as np                
from raysect.core import Point2D, Vector2D

####################################################################
####################################################################

def FirstOrderDerivative(f1 = None, f2 = None, dl = None):
    """
    Computes first order finite difference. Depending on points 1 and 2
    provided, backward, forward, or centred is implicitly employed.

    f1 and f2 are the values attained by the fuctions and are a distance
    dl apart in space. f2 follows f1.

    If dl <= 0 (coincident points or order reversed), then 0.0 is returned.
    """

    if dl > 0:
        return (f2 - f1) / dl
    else:
        return 0.0

####################################################################
####################################################################

def SecondOrderDerivative(f1 = None, f2 = None, f3 = None, f4 = None, dl1 = None, dl2 = None):
    """
    Computes second order finite difference. Depending on points 1, 2, 3 and 4
    provided, backward, forward, or centred is implicitly employed.

    f1 and f2 are the values attained by the fuctions and are a distance
    dl1 apart in space. f2 follows f1.

    f3 and f4 are the values attained by the fuctions and are a distance
    dl2 apart in space. f4 follows f3.

    If centred finite difference is in place, then f2 == f4 and the usual
    [-1 +2 -1] computational molecule is recovered.

    If dl1 <= 0 or dl2 <=0 (coincident points or order reversed), then 0.0 is returned.

    CAUTION.
    
    The factor (1 / dl1 / dl2) is a proxy (to avoid unnecessary complications) for the
    exact factor which would depend on the half-distance between points.
    This is "good enough" for the original purpose of the present function (non-uniform step).
    If this functions is employed in other context, this factor would need to be revisited.
    """

    if dl1 > 0 and dl2 > 0:
        return ((f2 - f1) - (f4 - f3)) / dl1 / dl2
    else:
        return 0.0

####################################################################
####################################################################

def ComputeGradient(cfg = None, b2fgmtry = None, SOLPSsim = None, is_homemade = None):
    """
    Given the SOLPS-ITER B2 structured mesh, computes the gradient vector
    in each cell of the ny x nx mesh (CHERAB convention - i.e transposed SOLPS-ITER
    convention).

    Where the gradient is non-defined, 0.0 is assumed. This is a proper assumption
    as long as the generation of the non-uniform step map is concerned. Must be
    revisited for other applications.
    """

    if is_homemade is False:

        # plasma quantities

#        if cfg['baserun'] == "radiation_load":
        f = SOLPSsim.total_radiation
#        elif cfg['baserun'] == "h_alpha_camera":
#            f = SOLPSsim.halpha_total_radiation
        
        mesh = SOLPSsim.mesh

        # mesh quantities
        
        nx = mesh.nx
        ny = mesh.ny
        cr = mesh.cr
        cz = mesh.cz

        neighbix = mesh.neighbix
        neighbiy = mesh.neighbiy

        leftix = neighbix[0,:,:]
        leftiy = neighbiy[0,:,:]

        bottomix = neighbix[1,:,:]
        bottomiy = neighbiy[1,:,:]

        rightix = neighbix[2,:,:]
        rightiy = neighbiy[2,:,:]

        topix = neighbix[3,:,:]
        topiy = neighbiy[3,:,:]

        # type of scenario

        if b2fgmtry['leftcut'].shape == (2,):
            is_double_null = True
            ix_midcut = int(b2fgmtry['leftcut'].sum())
        else:
            is_double_null = False

    elif is_homemade is True:

        is_double_null = False

        f = SOLPSsim['halpha_total_emission'].T

        nx = b2fgmtry['nx']
        ny = b2fgmtry['ny']
        cr = b2fgmtry['cx'].T
        cz = b2fgmtry['cy'].T

        leftix = b2fgmtry['leftix'].T
        leftiy = b2fgmtry['leftiy'].T

        bottomix = b2fgmtry['bottomix'].T
        bottomiy = b2fgmtry['bottomiy'].T

        rightix = b2fgmtry['rightix'].T
        rightiy = b2fgmtry['rightiy'].T

        topix = b2fgmtry['topix'].T
        topiy = b2fgmtry['topiy'].T

    else: raise ValueError('is_homemade must be specified!')

    # first order derivatives along x and y

    dfdx = np.zeros((ny, nx))
    dfdy = np.zeros((ny, nx))

    for ix in range(nx):
        for iy in range(ny):

            if is_double_null and ix == ix_midcut:

                # ad hoc because of double null mesh half-way split
                # => backward finite difference

                ix1 = leftix[iy,ix]
                ix2 = ix

            elif is_double_null and ix == ix_midcut + 1:

                # ad hoc because of double null mesh half-way split
                # => forward finite difference

                ix1 = ix
                ix2 = rightix[iy,ix]

            elif ix > 0 and ix < nx-1:

                # derivative along x
                # internal nodes
                # centred finite difference

                ix1 = leftix[iy,ix]
                ix2 = rightix[iy,ix]              

            elif ix == 0:

                # derivative along x
                # West nodes
                # forward finite difference

                ix1 = ix
                ix2 = rightix[iy,ix]

            elif ix == nx-1:

                # derivative along x
                # East nodes
                # backward finite difference

                ix1 = leftix[iy,ix]
                ix2 = ix

            f1 = f[iy,ix1]
            f2 = f[iy,ix2]
            c1 = Point2D(cr[iy,ix1], cz[iy,ix1])
            c2 = Point2D(cr[iy,ix2], cz[iy,ix2])
            dx = c1.distance_to(c2)
            
            dfdx[iy,ix] = FirstOrderDerivative(f1 = f1, f2 = f2, dl = dx)

            #################################

            if iy > 0 and iy < ny-1:

                # derivative along y
                # internal nodes
                # centred finite difference

                iy1 = bottomiy[iy,ix]
                iy2 = topiy[iy,ix]

            elif iy == 0:

                # derivative along y
                # South nodes
                # forward finite difference

                iy1 = iy
                iy2 = topiy[iy,ix]

            elif iy == ny-1:

                # derivative along y
                # North nodes
                # backward finite difference 

                iy1 = bottomiy[iy,ix]
                iy2 = iy

            f1 = f[iy1,ix]
            f2 = f[iy2,ix]
            c1 = Point2D(cr[iy1,ix], cz[iy1,ix])
            c2 = Point2D(cr[iy2,ix], cz[iy2,ix])
            dy = c1.distance_to(c2)
            
            dfdy[iy,ix] = FirstOrderDerivative(f1 = f1, f2 = f2, dl = dy)

    dfdx[np.isnan(dfdx)] = 0.0
    dfdx[np.isinf(dfdx)] = 0.0
    dfdy[np.isnan(dfdy)] = 0.0
    dfdy[np.isinf(dfdy)] = 0.0

    return [dfdx, dfdy]

####################################################################
####################################################################

def ComputeHessian(cfg = None, b2fgmtry = None, SOLPSsim = None, is_homemade = None):
    """
    Given the SOLPS-ITER B2 structured mesh, computes the Hessian matrix of the
    second order derivatives in each cell of the ny x nx mesh (CHERAB convention
    - i.e transposed SOLPS-ITER convention).

    Where the Hessian is non-defined, 0.0 is assumed. This is a proper assumption
    as long as the generation of the non-uniform step map is concerned. Must be
    revisited for other applications.
    """

    if is_homemade is False:

        # plasma quantities

#        if cfg['baserun'] == "radiation_load":
        f = SOLPSsim.total_radiation
#        elif cfg['baserun'] == "h_alpha_camera":
#            f = SOLPSsim.halpha_total_radiation
        
        mesh = SOLPSsim.mesh

        # mesh quantities
        
        nx = mesh.nx
        ny = mesh.ny
        cr = mesh.cr
        cz = mesh.cz

        neighbix = mesh.neighbix
        neighbiy = mesh.neighbiy

        leftix = neighbix[0,:,:]
        leftiy = neighbiy[0,:,:]

        bottomix = neighbix[1,:,:]
        bottomiy = neighbiy[1,:,:]

        rightix = neighbix[2,:,:]
        rightiy = neighbiy[2,:,:]

        topix = neighbix[3,:,:]
        topiy = neighbiy[3,:,:]

        # type of scenario

        if b2fgmtry['leftcut'].shape == (2,):
            is_double_null = True
            ix_midcut = int(b2fgmtry['leftcut'].sum())
        else:
            is_double_null = False

    elif is_homemade is True:

        is_double_null = False

        f = SOLPSsim['halpha_total_emission'].T

        nx = b2fgmtry['nx']
        ny = b2fgmtry['ny']
        cr = b2fgmtry['cx'].T
        cz = b2fgmtry['cy'].T

        leftix = b2fgmtry['leftix'].T
        leftiy = b2fgmtry['leftiy'].T

        bottomix = b2fgmtry['bottomix'].T
        bottomiy = b2fgmtry['bottomiy'].T

        rightix = b2fgmtry['rightix'].T
        rightiy = b2fgmtry['rightiy'].T

        topix = b2fgmtry['topix'].T
        topiy = b2fgmtry['topiy'].T

    else: raise ValueError('is_homemade must be specified!')

    # second order derivatives
    #
    # Schwarz theorem supposed to hold: SOLPS map smooth enough
    # so that d^2f/dxdy == d^2f/dydx for the mixed derivatives

    d2fdx2  = np.zeros((ny, nx))
    d2fdy2  = np.zeros((ny, nx))
    d2fdxdy = np.zeros((ny, nx))

    for ix in range(nx):
        for iy in range(ny):

            if is_double_null and ix == ix_midcut:

                # ad hoc because of double null mesh half-way split
                # => backward finite difference

                ix1 = leftix[iy,ix]
                ix2 = ix
                ix3 = leftix[iy,ix-1]
                ix4 = leftix[iy,ix]

            elif is_double_null and ix == ix_midcut + 1:

                # ad hoc because of double null mesh half-way split
                # => forward finite difference

                ix1 = rightix[iy,ix]
                ix2 = rightix[iy,ix+1]
                ix3 = ix
                ix4 = rightix[iy,ix]

            elif ix > 0 and ix < nx-1:

                # derivative along x
                # internal nodes
                # centred finite difference

                ix1 = ix
                ix2 = rightix[iy,ix]
                ix3 = leftix[iy,ix]
                ix4 = ix

            elif ix == 0:

                # derivative along x
                # West nodes
                # forward finite difference

                ix1 = rightix[iy,ix]
                ix2 = rightix[iy,ix+1]
                ix3 = ix
                ix4 = rightix[iy,ix]

            elif ix == nx-1:

                # derivative along x
                # East nodes
                # backward finite difference

                ix1 = leftix[iy,ix]
                ix2 = ix
                ix3 = leftix[iy,ix-1]
                ix4 = leftix[iy,ix]

            f1 = f[iy,ix1]
            f2 = f[iy,ix2]
            f3 = f[iy,ix3]
            f4 = f[iy,ix4]
            c1 = Point2D(cr[iy,ix1], cz[iy,ix1])
            c2 = Point2D(cr[iy,ix2], cz[iy,ix2])
            c3 = Point2D(cr[iy,ix3], cz[iy,ix3])
            c4 = Point2D(cr[iy,ix4], cz[iy,ix4])
            dl1 = c1.distance_to(c2)
            dl2 = c3.distance_to(c4)
            
            d2fdx2[iy,ix] = SecondOrderDerivative(f1 = f1, f2 = f2,
                                                  f3 = f3, f4 = f4,
                                                  dl1 = dl1, dl2 = dl2)

            #################################

            if iy > 0 and iy < ny-1:

                # derivative along y
                # internal nodes
                # centred finite difference

                iy1 = iy
                iy2 = topiy[iy,ix]
                iy3 = bottomiy[iy,ix]
                iy4 = iy

            elif iy == 0:

                # derivative along y
                # West nodes
                # forward finite difference

                iy1 = topiy[iy,ix]
                iy2 = topiy[iy+1,ix]
                iy3 = iy
                iy4 = topiy[iy,ix]

            elif iy == ny-1:

                # derivative along y
                # East nodes
                # backward finite difference

                iy1 = bottomiy[iy,ix]
                iy2 = iy
                iy3 = bottomiy[iy-1,ix]
                iy4 = bottomiy[iy,ix]

            f1 = f[iy1,ix]
            f2 = f[iy2,ix]
            f3 = f[iy3,ix]
            f4 = f[iy4,ix]
            c1 = Point2D(cr[iy1,ix], cz[iy1,ix])
            c2 = Point2D(cr[iy2,ix], cz[iy2,ix])
            c3 = Point2D(cr[iy3,ix], cz[iy3,ix])
            c4 = Point2D(cr[iy4,ix], cz[iy4,ix])
            dl1 = c1.distance_to(c2)
            dl2 = c3.distance_to(c4)

            d2fdy2[iy,ix] = SecondOrderDerivative(f1 = f1, f2 = f2,
                                                  f3 = f3, f4 = f4,
                                                  dl1 = dl1, dl2 = dl2)

            #################################

            if is_double_null and ix == ix_midcut and iy > 0 and iy < ny-1:

                # ad hoc because of double null mesh half-way split
                # => backward finite difference

                ix1 = ix
                ix2 = ix
                ix3 = leftix[iy,ix]
                ix4 = leftix[iy,ix]

                iy1 = bottomiy[iy,ix]
                iy2 = topiy[iy,ix]
                iy3 = bottomiy[iy,ix]
                iy4 = topiy[iy,ix]

            elif is_double_null and ix == ix_midcut + 1 and iy > 0 and iy < ny-1:

                # ad hoc because of double null mesh half-way split
                # => forward finite difference

                ix1 = rightix[iy,ix]
                ix2 = rightix[iy,ix]
                ix3 = ix
                ix4 = ix

                iy1 = bottomiy[iy,ix]
                iy2 = topiy[iy,ix]
                iy3 = bottomiy[iy,ix]
                iy4 = topiy[iy,ix]

            elif is_double_null and ix == ix_midcut and iy == 0:

                # ad hoc because of double null mesh half-way split
                # => backward finite difference

                ix1 = ix
                ix2 = ix
                ix3 = leftix[iy,ix]
                ix4 = leftix[iy,ix]

                iy1 = iy
                iy2 = topiy[iy,ix]
                iy3 = iy
                iy4 = topiy[iy,ix]

            elif is_double_null and ix == ix_midcut + 1 and iy == 0:

                # ad hoc because of double null mesh half-way split
                # => backward finite difference

                ix1 = rightix[iy,ix]
                ix2 = rightix[iy,ix]
                ix3 = ix
                ix4 = ix

                iy1 = iy
                iy2 = topiy[iy,ix]
                iy3 = iy
                iy4 = topiy[iy,ix]

            elif is_double_null and ix == ix_midcut and iy == ny-1:

                # ad hoc because of double null mesh half-way split
                # => forward finite difference

                ix1 = ix
                ix2 = ix
                ix3 = leftix[iy,ix]
                ix4 = leftix[iy,ix]

                iy1 = bottomiy[iy,ix]
                iy2 = iy
                iy3 = bottomiy[iy,ix]
                iy4 = iy

            elif is_double_null and ix == ix_midcut + 1 and iy == ny-1:

                # ad hoc because of double null mesh half-way split
                # => forward finite difference

                ix1 = rightix[iy,ix]
                ix2 = rightix[iy,ix]
                ix3 = ix
                ix4 = ix

                iy1 = bottomiy[iy,ix]
                iy2 = iy
                iy3 = bottomiy[iy,ix]
                iy4 = iy

            elif ix > 0 and ix < nx-1 and iy > 0 and iy < ny-1:

                # derivative along xy
                # internal nodes
                # centred finite difference

                ix1 = rightix[iy,ix]
                ix2 = rightix[iy,ix]
                ix3 = leftix[iy,ix]
                ix4 = leftix[iy,ix]

                iy1 = bottomiy[iy,ix]
                iy2 = topiy[iy,ix]
                iy3 = bottomiy[iy,ix]
                iy4 = topiy[iy,ix]

            elif ix == 0 and ix < nx-1 and iy > 0 and iy < ny-1:

                # derivative along xy
                # internal nodes
                # centred finite difference

                ix1 = rightix[iy,ix]
                ix2 = rightix[iy,ix]
                ix3 = ix
                ix4 = ix

                iy1 = bottomiy[iy,ix]
                iy2 = topiy[iy,ix]
                iy3 = bottomiy[iy,ix]
                iy4 = topiy[iy,ix]

            elif ix > 0 and ix == nx-1 and iy > 0 and iy < ny-1:

                # derivative along xy
                # internal nodes
                # centred finite difference

                ix1 = ix
                ix2 = ix
                ix3 = leftix[iy,ix]
                ix4 = leftix[iy,ix]

                iy1 = bottomiy[iy,ix]
                iy2 = topiy[iy,ix]
                iy3 = bottomiy[iy,ix]
                iy4 = topiy[iy,ix]

            elif ix > 0 and ix < nx-1 and iy == 0 and iy < ny-1:

                # derivative along xy
                # internal nodes
                # centred finite difference

                ix1 = rightix[iy,ix]
                ix2 = rightix[iy,ix]
                ix3 = leftix[iy,ix]
                ix4 = leftix[iy,ix]

                iy1 = iy
                iy2 = topiy[iy,ix]
                iy3 = iy
                iy4 = topiy[iy,ix]

            elif ix > 0 and ix < nx-1 and iy > 0 and iy == ny-1:

                # derivative along xy
                # internal nodes
                # centred finite difference

                ix1 = rightix[iy,ix]
                ix2 = rightix[iy,ix]
                ix3 = leftix[iy,ix]
                ix4 = leftix[iy,ix]

                iy1 = bottomiy[iy,ix]
                iy2 = iy
                iy3 = bottomiy[iy,ix]
                iy4 = iy

            elif ix == 0 and ix < nx-1 and iy == 0 and iy < ny-1:

                # derivative along xy
                # internal nodes
                # centred finite difference

                ix1 = rightix[iy,ix]
                ix2 = rightix[iy,ix]
                ix3 = ix
                ix4 = ix

                iy1 = iy
                iy2 = topiy[iy,ix]
                iy3 = iy
                iy4 = topiy[iy,ix]

            elif ix == 0 and ix < nx-1 and iy > 0 and iy == ny-1:

                # derivative along xy
                # internal nodes
                # centred finite difference

                ix1 = rightix[iy,ix]
                ix2 = rightix[iy,ix]
                ix3 = ix
                ix4 = ix

                iy1 = bottomiy[iy,ix]
                iy2 = iy
                iy3 = bottomiy[iy,ix]
                iy4 = iy

            elif ix > 0 and ix == nx-1 and iy == 0 and iy < ny-1:

                # derivative along xy
                # internal nodes
                # centred finite difference

                ix1 = ix
                ix2 = ix
                ix3 = leftix[iy,ix]
                ix4 = leftix[iy,ix]

                iy1 = iy
                iy2 = topiy[iy,ix]
                iy3 = iy
                iy4 = topiy[iy,ix]

            if ix > 0 and ix == nx-1 and iy > 0 and iy == ny-1:

                # derivative along xy
                # internal nodes
                # centred finite difference

                ix1 = ix
                ix2 = ix
                ix3 = leftix[iy,ix]
                ix4 = leftix[iy,ix]

                iy1 = bottomiy[iy,ix]
                iy2 = iy
                iy3 = bottomiy[iy,ix]
                iy4 = iy

            f1 = f[iy1,ix1]
            f2 = f[iy2,ix2]
            f3 = f[iy3,ix3]
            f4 = f[iy4,ix4]
            c1 = Point2D(cr[iy1,ix], cz[iy1,ix]) # CAUTION. Indeces...
            c2 = Point2D(cr[iy2,ix], cz[iy2,ix])
            c3 = Point2D(cr[iy3,ix3], cz[iy3,ix3])
            c4 = Point2D(cr[iy4,ix4], cz[iy4,ix4])
            dl1 = c1.distance_to(c2)
            dl2 = c3.distance_to(c4)
            
            d2fdxdy[iy,ix] = SecondOrderDerivative(f1 = f1, f2 = f2,
                                                   f3 = f3, f4 = f4,
                                                   dl1 = dl1, dl2 = dl2)

    d2fdx2[np.isnan(d2fdx2)] = 0.0
    d2fdx2[np.isinf(d2fdx2)] = 0.0
    d2fdy2[np.isnan(d2fdy2)] = 0.0
    d2fdy2[np.isinf(d2fdy2)] = 0.0
    d2fdxdy[np.isnan(d2fdxdy)] = 0.0
    d2fdxdy[np.isinf(d2fdxdy)] = 0.0

    return [d2fdx2, d2fdy2, d2fdxdy]

####################################################################
####################################################################

def ComputeGradientNorm(cfg = None, b2fgmtry = None, SOLPSsim = None, is_homemade = None):
    """
    Given the maps dfdx and dfdy of the first order derivatives,
    computes the norm of the gradient in each cell.
    """

    [dfdx, dfdy] = ComputeGradient(cfg = cfg, b2fgmtry = b2fgmtry, SOLPSsim = SOLPSsim, is_homemade = is_homemade)

    gradNorm = np.zeros(dfdx.shape)

    for ix in range(dfdx.shape[1]):
        for iy in range(dfdx.shape[0]):

            grad = Vector2D(dfdx[iy,ix], dfdy[iy,ix])
            gradNorm[iy,ix] = grad.length

    return gradNorm

####################################################################
####################################################################

def ComputeHessianNorm2(cfg = None, b2fgmtry = None, SOLPSsim = None, is_homemade = None):
    """
    Given the maps df2dx2, d2fdy2, and d2fdxdy of the second order derivatives,
    computes the 2-norm of the Hessian matrix in each cell (via eigenvalues).
    """

    [d2fdx2, d2fdy2, d2fdxdy] = ComputeHessian(cfg = cfg, b2fgmtry = b2fgmtry, SOLPSsim = SOLPSsim, is_homemade = is_homemade)

    hessianNorm2 = np.zeros(d2fdx2.shape)

    for ix in range(d2fdx2.shape[1]):
        for iy in range(d2fdx2.shape[0]):

            hessian = np.array([[d2fdx2[iy,ix],  d2fdxdy[iy,ix]],
                               [d2fdxdy[iy,ix], d2fdy2[iy,ix]]])
            matrix = np.matmul(np.transpose(hessian), hessian)
            hessianNorm2[iy,ix] = np.sqrt(np.linalg.eigvals(matrix).max())

    return hessianNorm2

####################################################################
####################################################################
