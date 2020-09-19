#!/usr/bin/env python

import math
import matplotlib.pyplot
import numpy

import libairfoil.parsec


def main():
       
    # jfoil_string = 'Parsec-11 [0.01:0.4:0.075:-0.1:0.4:-0.075:0.1:0:0:0:20]'
    # jfoil_string = 'Parsec-11 [0.06:0.3:0.118:-0.9:0.3:-0.118:0.9:0:0:0:20]'
    jfoil_string = 'Parsec-11 [0.0083:0.423:0.0587:-0.347:0.358:-0.032:0.417:0:0:10.03:5.64]'

    params = libairfoil.parsec.Parameters()
    params.load_from_javafoil_parsec11(jfoil_string)
    
    airfoil = libairfoil.parsec.Airfoil(params)
    
    matplotlib.pyplot.figure()
    x = numpy.linspace(0.0, 1.0, 200)
    matplotlib.pyplot.plot(x, airfoil.Z_up(x), 'r--', x, airfoil.Z_lo(x), 'b--')
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.show()
    

if __name__=='__main__':
    main()
