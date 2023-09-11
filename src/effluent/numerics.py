"""
The module contains numerical helper functions.
"""

import numpy as np


def bilin_inv(f, g, F, G, maxiter=7, tol=1.0e-7):
    """
    Inverse bilinear interpolation

    ``f, g`` should be scalars or arrays of same shape

    ``F, G`` should be 2D arrays of the same shape

    :param f: Desired f value
    :param g: Desired g value
    :param F: Tabulated f values
    :param G: Tabulated g values
    :param maxiter: Maximum number of Newton iterations
    :param tol: Maximum residual value
    :return: A tuple ``(x, y)`` such that ``F[x, y] = f`` and ``G[x, y] = g``, when
        linearly interpolated
    """

    imax, jmax = np.array(F.shape) - 1

    f = np.asarray(f)
    g = np.asarray(g)

    # initial guess
    x = np.zeros_like(f) + 0.5 * imax
    y = np.zeros_like(f) + 0.5 * jmax

    for t in range(maxiter):
        i = np.minimum(imax - 1, x.astype("i4"))
        j = np.minimum(jmax - 1, y.astype("i4"))

        p, q = x - i, y - j

        # Shorthands
        F00 = F[i, j]
        F01 = F[i, j+1]
        F10 = F[i+1, j]
        F11 = F[i+1, j+1]
        G00 = G[i, j]
        G01 = G[i, j+1]
        G10 = G[i+1, j]
        G11 = G[i+1, j+1]

        # Bilinear estimate of F[x,y] and G[x,y]
        Fs = (
            (1 - p) * (1 - q) * F00
            + p * (1 - q) * F10
            + (1 - p) * q * F01
            + p * q * F11
        )
        Gs = (
            (1 - p) * (1 - q) * G00
            + p * (1 - q) * G10
            + (1 - p) * q * G01
            + p * q * G11
        )

        H = (Fs - f) ** 2 + (Gs - g) ** 2

        if np.all(H < tol**2):
            break

        # Estimate Jacobi matrix
        Fx = (1 - q) * (F10 - F00) + q * (F11 - F01)
        Fy = (1 - p) * (F01 - F00) + p * (F11 - F10)
        Gx = (1 - q) * (G10 - G00) + q * (G11 - G01)
        Gy = (1 - p) * (G01 - G00) + p * (G11 - G10)

        # Newton-Raphson step
        # Jinv = np.linalg.inv([[Fx, Fy], [Gx, Gy]])
        # incr = - np.dot(Jinv, [Fs-f, Gs-g])
        # x = x + incr[0], y = y + incr[1]
        det = Fx * Gy - Fy * Gx
        x -= (Gy * (Fs - f) - Fy * (Gs - g)) / det
        y -= (-Gx * (Fs - f) + Fx * (Gs - g)) / det

        x = np.maximum(0, np.minimum(imax, x))
        y = np.maximum(0, np.minimum(jmax, y))

    return x, y
