"""
Contains the implementation for the Equation of State
"""

import numpy as np


def roms_rho(temp, salt, depth):
    """
    Computes water density from temperature, salinity and depth.

    The algorithm is taken directly from the ROMS source file ``rho_eos.F``

    David R. Jackett and Trevor J. Mcdougall (1995): |jackett1995|_.
    Journal of Atmospheric and Oceanic Technology 12, no. 2: 381–89.

    .. |jackett1995| replace:: Minimal Adjustment of Hydrographic Profiles to Achieve Static Stability
    .. _jackett1995: https://doi.org/10.1175/1520-0426(1995)012<0381:MAOHPT>2.0.CO;2

    :param temp: Temperature, in degrees Celcius
    :param salt: Salinity, in PSU
    :param depth: Depth, in meters
    :return: Density, in kg/m3
    """

    # Coefficients
    A00 = +1.909256e+04
    A01 = +2.098925e+02
    A02 = -3.041638e+00
    A03 = -1.852732e-03
    A04 = -1.361629e-05
    B00 = +1.044077e+02
    B01 = -6.500517e+00
    B02 = +1.553190e-01
    B03 = +2.326469e-04
    D00 = -5.587545e+00
    D01 = +7.390729e-01
    D02 = -1.909078e-02
    E00 = +4.721788e-01
    E01 = +1.028859e-02
    E02 = -2.512549e-04
    E03 = -5.939910e-07
    F00 = -1.571896e-02
    F01 = -2.598241e-04
    F02 = +7.267926e-06
    G00 = +2.042967e-03
    G01 = +1.045941e-05
    G02 = -5.782165e-10
    G03 = +1.296821e-07
    H00 = -2.595994e-07
    H01 = -1.248266e-09
    H02 = -3.508914e-09
    Q00 = +9.99842594e+02
    Q01 = +6.793952e-02
    Q02 = -9.095290e-03
    Q03 = +1.001685e-04
    Q04 = -1.120083e-06
    Q05 = +6.536332e-09
    U00 = +8.24493e-01
    U01 = -4.08990e-03
    U02 = +7.64380e-05
    U03 = -8.24670e-07
    U04 = +5.38750e-09
    V00 = -5.72466e-03
    V01 = +1.02270e-04
    V02 = -1.65460e-06
    W00 = +4.8314e-04

    #
    # Check temperature and salinity lower values. Assign depth to the pressure.
    #
    Tt = np.maximum(-2.5, temp)
    Tt = np.minimum(40, Tt)
    Ts = np.maximum(0, salt)
    Ts = np.minimum(100, Ts)
    sqrtTs = np.sqrt(Ts)

    Tp = -depth
    Tpr10 = 0.1 * Tp

    #
    # ----------------------------------------------------------------------------------
    # Compute density (kg/m3) at standard one atmosphere pressure.
    # ----------------------------------------------------------------------------------
    #
    C0 = Q00 + Tt * (Q01 + Tt * (Q02 + Tt * (Q03 + Tt * (Q04 + Tt * Q05))))
    C1 = U00 + Tt * (U01 + Tt * (U02 + Tt * (U03 + Tt * U04)))
    C2 = V00 + Tt * (V01 + Tt * V02)

    den1 = C0 + Ts * (C1 + sqrtTs * C2 + Ts * W00)

    #
    # ----------------------------------------------------------------------------------
    # Compute secant bulk modulus.
    # ----------------------------------------------------------------------------------
    #
    C3 = A00 + Tt * (A01 + Tt * (A02 + Tt * (A03 + Tt * A04)))
    C4 = B00 + Tt * (B01 + Tt * (B02 + Tt * B03))
    C5 = D00 + Tt * (D01 + Tt * D02)
    C6 = E00 + Tt * (E01 + Tt * (E02 + Tt * E03))
    C7 = F00 + Tt * (F01 + Tt * F02)
    C8 = G01 + Tt * (G02 + Tt * G03)
    C9 = H00 + Tt * (H01 + Tt * H02)

    bulk0 = C3 + Ts * (C4 + sqrtTs * C5)
    bulk1 = C6 + Ts * (C7 + sqrtTs * G00)
    bulk2 = C8 + Ts * C9
    bulk = bulk0 - Tp * (bulk1 - Tp * bulk2)

    #
    # ----------------------------------------------------------------------------------
    # Compute local "in situ" density anomaly (kg/m3 - 1000).
    # ----------------------------------------------------------------------------------
    #

    cff = 1.0 / (bulk + Tpr10)
    den = den1 * bulk * cff

    return den
