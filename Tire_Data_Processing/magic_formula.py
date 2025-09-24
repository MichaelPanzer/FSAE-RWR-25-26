"""
Contains functions for various versions of the Magic Formula tire model.
"""
import numpy as np

coordinate_mapping = {'ISO': (1,1,1,1,1,1)}

def mf61_Fy(x, p, scaling_factors=None, coordinate_system = 'adapted SAE'):
    """
    Pure lateral MF6.1 originally developed by TNO
    The equations below come from Tire and Vehicle Dynamics by Hans Pacejka 
    """
    #Model Parameters
    p_cy1, p_dy1, p_dy2, p_dy3, p_ey1, p_ey2, p_ey3, p_ey4, p_ey5, p_ky1, p_ky2, p_ky3, p_ky4, p_ky5, p_ky6, p_ky7, p_hy1, p_hy2, p_vy1, p_vy2, p_vy3, p_vy4, p_py1, p_py2, p_py4, p_py5, p_i0, f_z0 = p

    #Scaling Factors
    if scaling_factors is not None:
        lam_fz0, lam_cy, lam_muy, lam_ey, lam_kyalpha, lam_hy, lam_vy, lam_kygamma = scaling_factors
    else:
        lam_fz0 = lam_cy = lam_muy = lam_ey = lam_kyalpha = lam_hy = lam_vy = lam_kygamma = 1
    
    epsilon = 0.00


    #TODO ADD coordinate system transformations

        
    #Convert Model Inputs 
    alpha, F_z, gamma, p_i = x
    alpha = -np.radians(alpha)#NOTE this - sign is to convert the original 6.1 formula from adapted SAE to SAE coordinates
    F_z = np.abs(F_z)
    f_z0 = np.abs(f_z0)
    gamma = np.radians(gamma)

    #gamma_y = gamma * lam_gay  
    fz0pr  = f_z0 * lam_fz0 # 4.E1

    dfz = (F_z - fz0pr) / fz0pr  # 4.E2a
    dpi = (p_i-p_i0)/p_i0 # 4.E2b

    gamma_star = np.sin(gamma) # E4.4

    C_y = np.abs(p_cy1 * lam_cy)  # 4.E21 (>0) NOTE the abs might not be necessary but I think it is

    mu_y = (p_dy1 + p_dy2*dfz) * (1 + p_dy3*dpi + p_py4*dpi**2) * (1 - p_dy3*gamma_star**2) * lam_muy  # 4.E23
    D_y = mu_y * F_z  # 4.E22
    
    K_yalpha = p_ky1* fz0pr* (1 + p_py1*dpi)* (1-p_ky3*np.abs(gamma_star)) * np.sin( p_ky4* np.arctan(F_z / (fz0pr*(p_ky2 + p_ky5*gamma_star**2)*(1 + p_py2*dpi))) )  * lam_kyalpha  # 4.E25
    B_y = K_yalpha / (C_y * D_y)  # 4.E26 NOTE there might need to be a +epsilon_y in the denominator but I think that can be assumed to be 0

    S_vygamma = F_z * (p_vy3 + p_vy4*dfz) * gamma_star * lam_kygamma*lam_muy # E4.28
    S_vy = F_z*(p_vy1 + p_vy2*dfz)*lam_vy*lam_muy  + S_vygamma# 4.E29
    K_ygamma0 = F_z*(p_ky6 + p_ky7*dfz) * (1 + p_py5*dpi) * lam_kygamma# 4.E30

    S_hy = (p_hy1 + p_hy2*dfz) * lam_hy + (K_ygamma0*gamma_star-S_vygamma)/K_yalpha  # 4.E27
    alpha_y = alpha + S_hy  # 4.E20

    E_y = (p_ey1 + p_ey2*dfz) * (1 + p_ey5*gamma_star**2 - (p_ey3 + p_ey4*gamma_star)*np.sign(alpha_y)) * lam_ey  # 4.E24 (<=1)
    #E_y = [np.max(1.0,ey) for ey in E_y] #NOTE vectorized max function ill try this if it still doesn't work

    
    return D_y * np.sin(C_y * np.arctan(B_y*alpha_y - E_y*(B_y*alpha_y - np.arctan(B_y*alpha_y)))) + S_vy # 4.E19


def mf52_Fy(x, p, scaling_factors=None):
    # Example: Add a coordinate_system switch for future extension
    coordinate_system = 'default'  # or set this as a function argument if needed
    match coordinate_system:
        case 'default':
            pass  # Use default behavior
        case 'alternative':
            # Example: apply alternative coordinate system logic here
            pass
        case _:
            raise ValueError(f"Unknown coordinate system: {coordinate_system}")
    """
    Pure lateral MF5.2 originally developed by TNO
    This function has been adapted from Bill Cobb's matlab functions
    https://fsaettc.org/viewtopic.php?t=117
    """
    #Model Parameters
    p_cy1, p_dy1, p_dy2, p_dy3, p_ey1, p_ey2, p_ey3, p_ey4, p_ky1, p_ky2, p_ky3, p_hy1, p_hy2, p_hy3, p_vy1, p_vy2, p_vy3, p_vy4, fz0 = p

    #Scaling Factors
    if scaling_factors is not None:
        lam_fz0, lam_cx, lam_mux, lam_ex, lam_kx, lam_hx, lam_vx, lam_cy, lam_muy, lam_ey, lam_kyalpha, lam_hy, lam_vy, lam_gay, lam_tr, lam_res, lam_gaz, lam_xal, lam_yka, lam_vyka, lam_s, lam_sgkp, lam_sgal, lam_gyr = scaling_factors
    else:
        lam_fz0 =lam_cx = lam_mux = lam_ex = lam_kx = lam_hx = lam_vx = lam_cy = lam_muy = lam_ey = lam_kyalpha = lam_hy = lam_vy = lam_gay = lam_tr = lam_res = lam_gaz = lam_xal = lam_yka = lam_vyka = lam_s = lam_sgkp = lam_sgal = lam_gyr = 1


    #Convert Model Inputs 
    alpha, F_z, gamma = x
    alpha = -np.radians(alpha)
    F_z = np.abs(F_z)
    gamma = np.radians(gamma)

    gamma_y = gamma * lam_gay  # 31
    fz0pr  = fz0 * lam_fz0  # 15
    dfz    = (F_z - fz0pr) / fz0pr  # 14

    #-- lateral force (pure side slip)
    S_hy = (p_hy1+p_hy2 * dfz) * lam_hy + p_hy3 * gamma_y  # 38 lam_hy is new 
    alpha_y = alpha + S_hy  # 30
    C_y = p_cy1 * lam_cy  # 32
    mu_y = (p_dy1 + p_dy2 * dfz) * (1 - p_dy3*gamma_y**2) * lam_muy  # 34
    D_y = mu_y * F_z  # 33
    K_yalpha = p_ky1 * fz0 * np.sin(2 * np.arctan(F_z / (p_ky2 * fz0 * lam_fz0))) * (1-p_ky3 * np.abs(gamma_y)) * lam_fz0 * lam_kyalpha  # 36
    B_y = K_yalpha / (C_y * D_y)  # 37
    E_y = (p_ey1 + p_ey2 * dfz) * (1-(p_ey3 + p_ey4*gamma_y) * np.sign(alpha_y)) * lam_ey  # 35
    S_vy = F_z * ((p_vy1 + p_vy2*dfz)*lam_vy + (p_vy3 + p_vy4 * dfz) * gamma_y) * lam_muy  # 39 this is different from pacejka
    
    return D_y * np.sin(C_y * np.arctan(B_y * alpha_y - E_y*(B_y * alpha_y - np.arctan(B_y * alpha_y)))) + S_vy

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


