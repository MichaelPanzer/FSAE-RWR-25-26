"""
Contains functions for various versions of the Magic Formula tire model.
"""
import numpy as np

def mf_52():
    return None



def pacejka_lite_lateral(x, p):
    """
    5 term Pacejka-lite lateral tire model by Bill Cobb
    https://fsaettc.org/viewtopic.php?p=1777&hilit=pacejka+lite#p1777

    x: arraylike of inputs [alpha, Fz]
    p: arraylike of parameters [d1, d2, b, c, bp, SV]"""
    alpha, F_z = x

    if len(p)==6:
        D_1, D_2, B, C, B_p, S_v = p
    else:
        D_1, D_2, B, C, B_p = p
        S_v = 0

    F_z = -np.abs(F_z)  # Fz must always be negative!
    d = (D_1 + D_2 / 1000 * F_z) * F_z  # peak value (normalized)

    return d * np.sin((C + (B_p / 1000 * F_z)) * np.arctan(B * alpha)) + S_v


