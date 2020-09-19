import pytest

import numpy, math
import libairfoil.parsec as parsec


# def parse_javafoil_parsec11_params(string):
#     '''parses the JavaFoil parameter PARSEC-11 string into our PARSEC param object'''
#     jfoil = numpy.fromstring(string.replace(',', '.'), dtype=float, sep=':')
       
#     print(jfoil)
#     params = parsec.Parameters()
#     params.r_le       = jfoil[0]

#     params.X_lo       = jfoil[1]
#     params.Z_lo       = jfoil[2]
#     params.Z_XX_lo    = jfoil[3]
# What is this strings sequence ???
#     params.X_up       = jfoil[6]
#     params.Z_up       = jfoil[5]
#     params.Z_XX_up    = jfoil[4]

#     params.Z_te       = 0.0
#     params.dZ_te      = 0.0
#     params.alpha_te   = math.radians(0.0)
#     params.beta_te    = math.radians(20.0)

#     return params

def test_jfoil_sample():
    # parse jfoil string into our parameters
    # params = parse_javafoil_parsec11_params('0,01:0,4:0,075:-0,1: 0,075:0,075:-0,1:0:-0: 0:20')

    # todo: use code above after string parsing
    params = parsec.Parameters()
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

    airfoil = parsec.Airfoil(params)

    x_up = numpy.array([[1.00000000, 0.98116667, 0.95822222, 0.93150000, 0.90133333, 0.86805556,
                    0.83200000, 0.79350000, 0.75288889, 0.71050000, 0.66666667, 0.62172222,
                    0.57600000, 0.52983333, 0.48355556, 0.43750000, 0.39200000, 0.34738889,
                    0.30400000, 0.26216667, 0.22222222, 0.18450000, 0.14933333, 0.11705556,
                    0.08800000, 0.06250000, 0.04088889, 0.02350000, 0.01066667, 0.00272222,
                    0.00000000]])
    x_lo = numpy.array([[0.00000000, 0.00272222, 0.01066667, 0.02350000, 0.04088889, 0.06250000,
                    0.08800000, 0.11705556, 0.14933333, 0.18450000, 0.22222222, 0.26216667,
                    0.30400000, 0.34738889, 0.39200000, 0.43750000, 0.48355556, 0.52983333,
                    0.57600000, 0.62172222, 0.66666667, 0.71050000, 0.75288889, 0.79350000,
                    0.83200000, 0.86805556, 0.90133333, 0.93150000, 0.95822222, 0.98116667,
                    1.00000000]])
    
    z_up_expected = numpy.array([[0.00000000, 0.00365591, 0.00883517, 0.01550735, 0.02341024, 0.03208125,
                    0.04093427, 0.04936126, 0.05683477, 0.06298945, 0.06766698, 0.07091725, 0.07295802,
                    0.07410356, 0.07467861, 0.07493622, 0.07499670, 0.07482051, 0.07422028, 0.07290934,
                    0.07057594, 0.06696678, 0.06196035, 0.05561203, 0.04815822, 0.03997501, 0.03149790,
                    0.02311931, 0.01508883, 0.00744447, 0.00000000]])
    
    z_lo_expected = numpy.array([[0.00000000, -0.00744447, -0.01508883, -0.02311931, -0.03149790, -0.03997501,
                    -0.04815822, -0.05561203, -0.06196035, -0.06696678, -0.07057594, -0.07290934, -0.07422028,
                    -0.07482051, -0.07499670, -0.07493622, -0.07467861, -0.07410356, -0.07295802, -0.07091725,
                    -0.06766698, -0.06298945, -0.05683477, -0.04936126, -0.04093427, -0.03208125, -0.02341024,
                    -0.01550735, -0.00883517, -0.00365591, 0.00000000]])
   
    z_up = airfoil.Z_up(x_up)
    assert numpy.allclose(z_up_expected, z_up)
