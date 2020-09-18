'''
#          Copyright Marc Bodmer 2013.
# Distributed under the Boost Software License, Version 1.0.
#    (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)


Python library for parametric representation of airfoils according to PARSEC (PARametric SECtion).

References:

- Parametric Airfoils and Wings, by Helmut Sobieczky,
  http://www.as.dlr.de/hs/h-pdf/H141.pdf
  
- Representation Method Effects on Vibrational Genetic Algorithm in 2D Airfoil Design, Y.V.Pehlivanoglu
  http://www.hho.edu.tr/HutenDergi/2009TEMMUZ/5_PEHLIVANOGLU.pdf

- Aerodynamic Shape Optimization using Computer Mapping of Natural Evolution Process, Selvakumar/Mukesh
  http://sumo.intec.ugent.be/system/files/private/Selva_1_china.pdf
'''

import numpy
import math

class Parameters(object):
    '''Parameters defining a PARSEC airfoil'''   
    def __init__(self):
        self.r_le       = 0.0   # Leading edge radius
        self.X_up       = 0.0   # Upper crest location X coordinate
        self.Z_up       = 0.0   # Upper crest location Z coordinate
        self.Z_XX_up    = 0.0   # Upper crest location curvature
        self.X_lo       = 0.0   # Lower crest location X coordinate
        self.Z_lo       = 0.0   # Lower crest location Z coordinate
        self.Z_XX_lo    = 0.0   # Lower crest location curvature
        self.Z_te       = 0.0   # Trailing edge Z coordinate
        self.dZ_te      = 0.0   # Trailing edge thickness
        self.alpha_te   = 0.0   # Trailing edge direction angle
        self.beta_te    = 0.0   # Trailing edge wedge angle
        self.P_mix      = 1.0   # Blending parameter


class Coefficients(object):
    '''
    This class calculates the equation systems which define the coefficients
    for the polynomials given by the parsec airfoil parameters.
    '''
    def __init__(self, parsec_params):
        self._a_up = self._calc_a_up(parsec_params)
        self._a_lo = self._calc_a_lo(parsec_params)
    
    def a_up(self):
        '''Returns coefficient vector for upper surface'''
        return self._a_up
    
    def a_lo(self):
        '''Returns coefficient vector for lower surface'''
        return self._a_lo
    
    def _calc_a_up(self, parsec_params):
        Amat = self._prepare_linsys_Amat(parsec_params.X_up)
        Bvec = numpy.array([parsec_params.Z_te, parsec_params.Z_up,
                            math.tan(parsec_params.alpha_te - parsec_params.beta_te/2),
                            0.0, parsec_params.Z_XX_up, math.sqrt(2*parsec_params.r_le)]) 
        return numpy.linalg.solve(Amat, Bvec)
    
    def _calc_a_lo(self, parsec_params):
        Amat = self._prepare_linsys_Amat(parsec_params.X_lo)
        Bvec = numpy.array([parsec_params.Z_te, parsec_params.Z_lo,
                            math.tan(parsec_params.alpha_te + parsec_params.beta_te/2),
                            0.0, parsec_params.Z_XX_lo, -math.sqrt(2*parsec_params.r_le)])
        return numpy.linalg.solve(Amat, Bvec)
    
    def _prepare_linsys_Amat(self, X):
        return numpy.array(
            [[1.0,           1.0,          1.0,         1.0,          1.0,          1.0        ],
             [X**0.5,        X**1.5,       X**2.5,      X**3.5,       X**4.5,       X**5.5     ],
             [0.5,           1.5,          2.5,         3.5,          4.5,          5.5        ],
             [0.5*X**-0.5,   1.5*X**0.5,   2.5*X**1.5,  3.5*X**2.5,   4.5*X**3.5,   5.5*X**4.5 ],
             [-0.25*X**-1.5, 0.75*X**-0.5, 3.75*X**0.5, 8.75*X**1.5, 15.75*X**2.5, 24.75*X**3.5],
             [1.0,           0.0,          0.0,         0.0,          0.0,          0.0        ]])


class Airfoil(object):
    '''Airfoil defined by PARSEC Parameters'''
    def __init__(self, parsec_params):
        self._coeff = Coefficients(parsec_params)
        
    def Z_up(self, X):
        '''Returns Z(X) on upper surface, calculates PARSEC polynomial'''
        a = self._coeff.a_up()
        print(a)
        return a[0]*X**0.5 + a[1]*X**1.5 + a[2]*X**2.5 + a[3]*X**3.5 + a[4]*X**4.5 + a[5]*X**5.5
        
    
    def Z_lo(self, X):
        '''Returns Z(X) on lower surface, calculates PARSEC polynomial'''
        a = self._coeff.a_lo()
        print(a)
        return a[0]*X**0.5 + a[1]*X**1.5 + a[2]*X**2.5 + a[3]*X**3.5 + a[4]*X**4.5 + a[5]*X**5.5
