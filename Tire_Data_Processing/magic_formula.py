"""
Contains functions for various versions of the Magic Formula tire model.
"""
import numpy as np

def pacejka_lite_lateral(x, D_1, D_2, B, C, B_p):
    """
    5 term Pacejka-lite lateral tire model by Bill Cobb
    https://fsaettc.org/viewtopic.php?p=1777&hilit=pacejka+lite#p1777

    p: arraylike of parameters [d1, d2, b, c, bp, SV]
    x: arraylike of inputs [alpha, Fz]"""
    alpha, F_z = x

    F_z = -np.abs(F_z)  # Fz must always be negative!
    d = (D_1 + D_2 / 1000 * F_z) * F_z  # peak value (normalized)

    return d * np.sin((C + (B_p / 1000 * F_z)) * np.arctan(B * alpha))


