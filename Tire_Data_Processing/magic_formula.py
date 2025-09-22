"""
Contains functions for various versions of the Magic Formula tire model.
"""
import numpy as np

def mf52_Fy(x, p, scaling_factors=None):
    """
    Pure lateral MF5.2 originally developed by TNO
    This function has been adapted from Bill Cobb's matlab functions
    https://fsaettc.org/viewtopic.php?t=117
    """
    #Model Parameters
    p_cy1, p_dy1, p_dy2, p_dy3, p_ey1, p_ey2, p_ey3, p_ey4, p_ky1, p_ky2, p_ky3, p_hy1, p_hy2, p_hy3, p_vy1, p_vy2, p_vy3, p_vy4, p_fz0 = p

    #Scaling Factors
    if scaling_factors is not None:
        lam_f0, lam_cx, lam_mux, lam_ex, lam_kx, lam_hx, lam_vx, lam_cy, lam_muy, lam_ey, lam_ky, lam_hy, lam_vy, lam_gay, lam_tr, lam_res, lam_gaz, lam_xal, lam_yka, lam_vyka, lam_s, lam_sgkp, lam_sgal, lam_gyr = scaling_factors
    else:
        lam_f0, lam_cx = lam_mux = lam_ex = lam_kx = lam_hx = lam_vx = lam_cy = lam_muy = lam_ey = lam_ky = lam_hy = lam_vy = lam_gay = lam_tr = lam_res = lam_gaz = lam_xal = lam_yka = lam_vyka = lam_s = lam_sgkp = lam_sgal = lam_gyr = 1


    #Model Inputs 
    alpha, fz, gamma = x

    alpha = np.radians(alpha)
    fz = np.abs(fz)
    gamma = np.radians(gamma)


    gamma_y = gamma * lam_gay  # 31
    fz0pr  = p_fz0 * lam_f0  # 15
    dfz    = (fz - fz0pr) / fz0pr  # 14

    #-- lateral force (pure side slip)
    shy = (p_hy1 + p_hy2 * dfz) * lam_hy + p_hy3 * gamma_y  # 38 lam_hy is new 
    alpha_y = alpha + shy  # 30
    cy = p_cy1 * lam_cy  # 32
    muy = (p_dy1 + p_dy2*dfz) * (1.0 - p_dy3*gamma_y**2) * lam_muy  # 34
    dy = muy * fz  # 33
    ky = p_ky1 * p_fz0 * np.sin(2.0 * np.arctan(fz / (p_ky2 * p_fz0 * lam_f0))) * (1.0 - p_ky3 * np.abs(gamma_y)) * lam_f0 * lam_ky  # 36
    by = ky / (cy * dy)  # 37
    ey = (p_ey1 + p_ey2 * dfz) * (1.0 - (p_ey3 + p_ey4 * gamma_y) * np.sign(alpha_y)) * lam_ey  # 35
    svy = fz * ((p_vy1 + p_vy2 * dfz) * lam_vy + (p_vy3 + p_vy4 * dfz) * gamma_y) * lam_muy  # 39
    
    return dy * np.sin(cy * np.arctan(by * alpha_y - ey * (by * alpha_y - np.arctan(by * alpha_y)))) + svy




def pacejka_lite_Fy(x, p):
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


