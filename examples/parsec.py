#!/usr/bin/env python

import math
import matplotlib.pyplot
import numpy

import libairfoil.parsec


def main():
       
    params = libairfoil.parsec.Parameters()
#     params.r_le       = 0.02
#     params.X_up       = 0.35
#     params.Z_up       = 0.08
#     params.Z_XX_up    = -0.65
#     params.X_lo       = 0.25
#     params.Z_lo       = -0.08
#     params.Z_XX_lo    = 0.6
#     params.Z_te       = 0.0
#     params.dZ_te      = 0.04
#     params.alpha_te   = math.radians(-5.0)
#     params.beta_te    = math.radians(8.0)

#R_LE;  x_PRE;y_PRE;d2y/dx2_PRE;th_PRE;  x_SUC;y_SUC;d2y/dx2_SUC;th_SUC
#0.01;  0.450;-0.006;-0.2;0.05;          0.350;0.055;-0.350;-6

#Parsec-11 [0,01: 0,4: 0,075: -0,1: 0,075: 0,075: -0,1: 0: -0: 0: 20]

    params.r_le       = 0.01
    params.X_up       = 0.4
    params.Z_up       = 0.075
    params.Z_XX_up    = -0.1
    params.X_lo       = 0.4
    params.Z_lo       = -0.075
    params.Z_XX_lo    = 0.1
    params.Z_te       = 0.0
    params.dZ_te      = 0.0
    params.alpha_te   = math.radians(0.0)
    params.beta_te    = math.radians(20.0)
    
    airfoil = libairfoil.parsec.Airfoil(params)
    
    matplotlib.pyplot.figure()
    x = numpy.linspace(0.0, 1.0, 200)
    matplotlib.pyplot.plot(x, airfoil.Z_up(x), 'r--', x, airfoil.Z_lo(x), 'b--')
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.show()
    

if __name__=='__main__':
    main()
