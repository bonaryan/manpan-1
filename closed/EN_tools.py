from math import sqrt, pi, log

###### SIMPLY SUPPORTED PLATE ######

def N_pl_Rd(thickness, width, f_yield, psi=None):
    # Convert inputs to floats
    thickness, width, f_yield = float(thickness), float(width), float(f_yield)
    
    # Elastic critical load acc. to EN3-1-5 Annex A
    # Uniform compression is assumed (k_sigma = 4)
    if psi is None:
        psi = 1.
    else:
        psi = float(psi)
    
    # Calculate kapa_sigma
    kapa_sigma = 8.2 / (1.05 + psi)
    
    # Critical compression load
    N_cr_plate = thickness * width * kapa_sigma * sigma_cr_plate(thickness, width)
	
    # Aeff calculation.
    # Reduction factor for the effective area of the profile acc. to EC3-1-5
    classification = width / (thickness * sqrt(235 / f_yield))
    lambda_p = classification / (28.4 * sqrt(kapa_sigma))
    if lambda_p > 0.673 and plate_class(width, thickness, f_yield) == 4:
        rho = (lambda_p - 0.055 * (3 + psi)) / lambda_p ** 2
    else:
        rho = 1.
    
    # Effective area
    A_eff = rho * thickness * width
    
    # Axial compression resistance , Npl
    N_pl_Rd = A_eff * f_yield
    
    # Return value
    return N_pl_Rd


def plate_class(width, thickness, f_yield):
    # Convert inputs to floats
    width, thickness, f_yield = float(width), float(thickness), float(f_yield)
    
    # Calculate classification
    classification = width / (thickness * sqrt(235. / f_yield))
    if  classification <= 33.:
        p_class = 1
    elif classification <= 38.:
        p_class = 2
    elif classification <= 42.:
        p_class = 3
    else:
        p_class = 4
    
    # Return value
    return p_class


def sigma_cr_plate(thickness, width):
    # Convert inputs to floats
    thickness, width = float(thickness), float(width)
    
    # Elastic critical stress acc. to EN3-1-5 Annex A
    sigma_cr =  190000 * (thickness / width) ** 2
    
    # Return value
    return sigma_cr


###### CYLINDRICAL SHELLS ######

def sigma_x_Rd(thickness, radius, length, f_y_k, fab_quality = None, gamma_M1 = None):
    # Default values
    if fab_quality is None:
        fab_quality = 'fcA'
    elif not((fab_quality is 'fcA') or (fab_quality is 'fcB') or (fab_quality is 'fcC')):
        print('Invalid fabrication class input. Choose between \'fcA\', \'fcB\' and \'fcC\' ')
    
    if gamma_M1 is None:
        gamma_M1 = 1.
    else:
        gamma_M1 = float(gamma_M1)
    
    # Fabrication quality class acc. to table D2
    if fab_quality is 'fcA':
        Q_factor = 40.
    elif fab_quality == 'fcB':
        Q_factor = 25.
    elif fab_quality == 'fcC':
        Q_factor = 16.
    
    # Critical meridinal stress, calculated on separate function
    sigma_cr = sigma_x_Rcr(thickness, radius, length)
    
    # Shell slenderness
    lmda = sqrt(f_y_k / sigma_cr[0])
    delta_w_k = (1. / Q_factor) * sqrt(radius / thickness) * thickness
    alpha = 0.62 / (1 + 1.91 * (delta_w_k / thickness) ** 1.44)
    beta = 0.6
    eta = 1.
    if sigma_cr[1] is 'long':
        # For long cylinders, a formula is suggested fo lambda, EC3-1-6 D1.2.2(4)
        # Currently, the general form is used. to be fixed.
        lmda_0 = 0.2
        #lmda_0 = 0.2 + 0.1 * (sigma_E_M / sigma_E)
    else:
        lmda_0 = 0.2
    
    lmda_p = sqrt(alpha / (1. - beta))
    
    # Buckling reduction factor, chi
    if lmda <= lmda_0:
        chi = 1.
    elif lmda < lmda_p:
        chi = 1. - beta * ((lmda - lmda_0) / (lmda_p - lmda_0)) ** eta
    else:
        chi = alpha / (lmda ** 2)
    
    # Buckling stress
    sigma_Rk = chi * f_y_k
    sigma_Rd = sigma_Rk / gamma_M1
    
    # Return value
    return sigma_Rd


def N_cr_shell(thickness, radius, length):
    # Convert inputs to floats
    thickness, radius, length = float(thickness), float(radius), float(length)
    
    # Elastic critical load acc to EN3-1-6 Annex D
    N_cr_shell = 2 * pi * radius * thickness * sigma_x_Rcr(thickness, radius, length)[0]
    
    # Return value
    return N_cr_shell


def sigma_x_Rcr(thickness, radius, length):
    # Convert inputs to floats
    thickness, radius, length = float(thickness), float(radius), float(length)
    
    # Elastic critical load acc. to EN3-1-6 Annex D
    omega = length / sqrt(radius * thickness)
    if 1.7 <= omega and omega <= 0.5 * (radius / thickness):
        C_x = 1.
        length_category = 'medium'
    elif omega < 1.7:
        C_x = 1.36 - (1.83 / omega) + (2.07 / omega ** 2)
        length_category = 'short'
    else:
        # C_x_b is read on table D.1 of EN3-1-5 Annex D acc. to BCs
        # BC1 - BC1 is used on the Abaqus models (both ends clamped, see EN3-1-5 table 5.1)
        C_x_b = 6.
        C_x_N = max((1 + 0.2 * (1 - 2 * omega * thickness / radius) / C_x_b), 0.6)
        C_x = C_x_N
        length_category = 'long'
    
    # Calculate critical stress, eq. D.2 on EN3-1-5 D.1.2.1-5
    sigma_cr = 0.605 * 210000 * C_x * thickness / radius
    
    # Return value
    return sigma_cr, length_category


###### FLEXURAL BUCKLING ######

def lmbda(
    length,
    area,
    I_2,
    kapa_BC = None,
    E_modulus = None,
    f_yield = None
    ):
    
    # default values
    if kapa_BC is None:
        kapa_BC = 1.
    else:
        kapa_BC = float(kapa_BC)
    
    if E_modulus is None:
        E_modulus = 210000.
    else:
        E_modulus = float(E_modulus)
    
    if f_yield is None:
        f_yield = 380.
    else:
        f_yield = float(f_yield)
    
    # Calculate Euler's critical load
    N_cr = N_cr_flex(
        length,
        I_2,
        E_modulus = E_modulus,
        kapa_BC = kapa_BC
        )
    
    # Flexural slenderness EN3-1-1 6.3.1.3 (1)
    lmbda = sqrt(area * f_yield / N_cr)
    
    # Return the result
    return lmbda


def N_cr_flex(
    length,
    I_2,
    E_modulus = None,
    kapa_BC = None
    ):
    
    # default values
    if kapa_BC is None:
        kapa_BC = 1.
    else:
        kapa_BC = float(kapa_BC)
    
    if E_modulus is None:
        E_modulus = 210000.
    else:
        E_modulus = float(E_modulus)
    
    # Euler's critical load
    N_cr_flex = (pi ** 2) * E_modulus * I_2 / (kapa_BC * length) ** 2
    
    # Return the result
    return N_cr_flex


###### CIRCULAR PLATE  #####
# The following formulas are taken from table 18-3 of:
# W. D. Pilkey, Formulas for stress, strain, and structural matrices. New York: Wiley, 1994.
#
## Hinged perimeter
### Concentrated force applied on circle
#############(not finished, there is an errror, to be revised)###########

def circular_plate(a_L, a_1, W_load, r):
    poisson = 0.3
    
    alfa = r / a_L
    beta = a_1 / a_L
    
    C_1 = (3 + poisson) / (1 + poisson) * (1 + beta ** 2) + 2 * beta ** 2 * log(beta)
    C_2 = (1 - poisson) / (1 + poisson) * (1 - beta ** 2) + 2 * log(beta)
    C_3 = (3 + poisson) / (1 + poisson) - beta ** 2 * (1 - poisson) / (1 + poisson)
    
    if alfa <= beta:
        M_r = (1 / 4) * W_load * a_L * (1 + poisson) * C_2
        M_phi = M_r
    else:
        M_r = (1 / 4) * (W_load * a_1) / (alfa ** 2) * ((1 - poisson) * (1 - alfa ** 2) * beta ** 2 - 2 * (1 + poisson) * alfa ** 2 * log(alfa))
        M_phi = (1 / 4) * (W_load * a_1) / (alfa ** 2) * (2 * (1 - poisson) * alfa ** 2 - (1 - poisson) * (1 + alfa ** 2) * beta ** 2 - 2 * (1 + poisson) * alfa ** 2 * log(alfa))
    
    return M_r, M_phi