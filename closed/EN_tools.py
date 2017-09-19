"""This is the general docstring"""
from math import sqrt, pi, log
###### SIMPLY SUPPORTED PLATE ######


def N_pl_Rd(
            thickness,
            width,
            f_yield,
            psi
            ):
    # Docstring
    """
    Plastic design resistance of a plate.
    
    Calculates the resistance of a plate according to EN1993-1-1 and 
    EN1993-1-5. The plate is assumed simply supported.
    
    Parameters
    ----------
    thickness : float
        [mm] Plate thickness
    width : float
        [mm] Plate width
    f_yield : float
        [MPa] Yield stress
    psi : float, optional
        [_] Ratio of the min over max stress for a linear distribution, 
        (sigma_min / sigma_max)
        Default = 1, which implies a uniform distribution
    
    Returns
    -------
    float
        [N] Plastic design resistance
    
    Notes
    -----
    To be extended to include cantilever plate (outstand members)
    
    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings. Brussels: CEN, 2005.
    .. [2] Eurocode 3: Design of steel structures - Part 1-5: Plated structural elements. Brussels: CEN, 2005.

    """
    
    # Convert inputs to floats
    thickness, width, f_yield = float(thickness), float(width), float(f_yield)
    
    # Default value for psi
    if psi is None:
        psi = 1.
    else:
        psi = float(psi)
    
    # Calculate kapa_sigma
    k_sigma = 8.2 / (1.05 + psi)

    # Critical compression load
    N_cr_plate = thickness * width * sigma_cr_plate(thickness, width, psi=psi)
    
    # Aeff calculation.
    # Reduction factor for the effective area of the profile acc. to EC3-1-5
    classification = width / (thickness * sqrt(235 / f_yield))
    lambda_p = classification / (28.4 * sqrt(k_sigma))
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


def plate_class(
                thickness,
                width,
                f_yield
                ):
    
    # Docstring
    """
    Plate classification.
    
    Returnes the class for a given plate, according to EN1993-1-1.
    Currently works for simply supported plates under pure compression.
    
    Parameters
    ----------
    thickness : float
        [mm] Plate thickness
    width : float
        [mm] Plate width
    f_yield : float
        [MPa] Yield stress
    
    Returns
    -------
    int
        [_] Class number
    
    Notes
    -----
    To be extended to include the rest of the cases of Table 5.3 [1].
    Members under combined axial and bending and outstand members.
    
    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings. Brussels: CEN, 2005.
    
    """
    

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


def sigma_cr_plate(
                   thickness, 
                   width, 
                   psi=None
                   ):
    
    # Docstring
    """
    Critical stress of a plate.
    
    Calculates the critical stress for a simply supported plate.
    
    Parameters
    ----------
    thickness : float
        [mm] Plate thickness
    width : float
        [mm] Plate width
    psi : float, optional
        [_] Ratio of the min over max stress for a linear distribution, 
        (sigma_min / sigma_max)
        Default = 1, which implies a uniform distribution
    
    Returns
    -------
    float
        [MPa] Plate critical stress
    
    Notes
    -----
    To be extended to include cantilever plate (outstand members)
    
    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-5: Plated structural elements. Brussels: CEN, 2005.

    """
    # Convert inputs to floats
    thickness, width = float(thickness), float(width)
    
    # Default value for psi
    if psi is None:
        psi = 1.
    else:
        psi = float(psi)
    
    # Calculate kapa_sigma
    k_sigma = 8.2 / (1.05 + psi)
    
    # Elastic critical stress acc. to EN3-1-5 Annex A
    sigma_E =  190000 * (thickness / width) ** 2
    sigma_cr = sigma_E * k_sigma
    
    # Return value
    return sigma_cr


###### CYLINDRICAL SHELLS ######

def sigma_x_Rd(
               thickness,
               radius,
               length,
               f_y_k,
               fab_quality = None,
               gamma_M1 = None
               ):
    
    # Docstring
    """ 
    Meridional design buckling stress.
    
    Calculates the meridional buckling stress for a cylindrical shell 
    according to EN1993-1-6 [1].
    
    Parameters
    ----------
    thickness : float
        [mm] Shell thickness
    radius : float
        [mm] Cylinder radius
    length : float
        [mm] Cylnder length
    f_y_k : float
        [MPa] Characteristic yield strength
    fab_quality : str, optional
        [_] Fabrication quality class. Accepts: 'fcA', 'fcB', 'fcC'
        The three classes correspond to .006, .010 and .016 times the 
        width of a dimple on the shell.
        Default = 'fcA', which implies excelent fabrication
    gamma_M1 : int, optional
        [_] Partial safety factor
        Default = 1.
    
    Returns
    -------
    float
        [MPa] Meridional buckling stress
    
    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-6: Strength and stability of shell structures. Brussels: CEN, 2006.

    """

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


def N_cr_shell(
               thickness, 
               radius, 
               length
               ):
    
    # Docstring
    """ 
    Critical compressive load for cylindrical shell.
    
    Calculates the critical load for a cylindrical shell under pure 
    compression and assumes uniform stress distribution. Calculation
    according to EN1993-1-6 [1], Annex D.
    
    Parameters
    ----------
    thickness : float
        [mm] Shell thickness
    radius : float
        [mm] Cylinder radius
    length : float
        [mm] Cylnder length
    
    Returns
    -------
    float
        [N] Critical load
    
    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-6: Strength and stability of shell structures. Brussels: CEN, 2006.

    """

    # Convert inputs to floats
    thickness, radius, length = float(thickness), float(radius), float(length)
    
    # Elastic critical load acc to EN3-1-6 Annex D
    N_cr_shell = 2 * pi * radius * thickness * sigma_x_Rcr(thickness, radius, length)[0]
    
    # Return value
    return N_cr_shell


def sigma_x_Rcr(
                thickness, 
                radius, 
                length
                ):
    
    # Docstring
    """ 
    Critical meridional stress for cylindrical shell.
    
    Calculates the critical load for a cylindrical shell under pure 
    compression and assumes uniform stress distribution. Calculation
    according to EN1993-1-6 [1], Annex D.
    
    Parameters
    ----------
    thickness : float
        [mm] Shell thickness
    radius : float
        [mm] Cylinder radius
    length : float
        [mm] Cylnder length
    
    Returns
    -------
    float
        [N] Critical load
    
    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-6: Strength and stability of shell structures. Brussels: CEN, 2006.

    """
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
    MOI,
    kapa_BC = None,
    E_modulus = None,
    f_yield = None
    ):
    
    # Docstring
    """
    Flexural slenderness.
    
    Calculates the slenderness of a columne under pure compression.
    Euler's critical load is used.
    
    Parameters
    ----------
    length : float
        [mm] Column length
    area : float
        [mm^2] Cross section area
    MOI : float
        [mm^4] Moment of inertia
    kapa_BC : float, optional
        [_] length correction for the effect of the boundary conditions.
        Default = 1, which implies simply supported column
    E_modulus : float, optional
        [MPa] Modulus of elasticity
        Default = 210000., typical value for steel
    f_yield : float, optional
        [MPa] yield stress.
        Default = 380., brcause this value was used extencively while the
        function was being written. To be changed to 235.
    
    Returns
    -------
    float
        [_] Member slenderness
    
    """
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
        MOI,
        E_modulus = E_modulus,
        kapa_BC = kapa_BC
        )
    
    # Flexural slenderness EN3-1-1 6.3.1.3 (1)
    lmbda = sqrt(area * f_yield / N_cr)
    
    # Return the result
    return lmbda


def N_cr_flex(
    length,
    MOI,
    kapa_BC = None,
    E_modulus = None
    ):
    
    # Docstring
    """
    Euler's critical load.
    
    Calculates the critical load for flexural buckling of a given column.
        
    Parameters
    ----------
    length : float
        [mm] Column length.
    MOI : float
        [mm^4] Moment of inertia.
    kapa_BC : float, optional
        [_] length correction for the effect of the boundary conditions.
        Default = 1, which implies simply supported column.
    E_modulus : float, optional
        [MPa] Modulus of elasticity.
        Default = 210000., typical value for steel.
    
    Returns
    -------
    float
        [N] Critical load.
    
    """
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
    N_cr_flex = (pi ** 2) * E_modulus * MOI / (kapa_BC * length) ** 2
    
    # Return the result
    return N_cr_flex


def imp_factor(b_curve):
    
    # Docstring
    """
    Imperfection factor.
    
    Returns the imperfection factor for a given buckling curve.
    The values are taken from Table 6.1 of EN1993-1-1 [1]
        
    Parameters
    ----------
    b_curve : str
        [_] Name of the buckling curve as obtained from Table 6.2 of [1].
        Valid options are {'a0', 'a', 'b', 'c', 'd'}
    
    Returns
    -------
    float
        [_] Imperfection factor.
    
    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings. Brussels: CEN, 2005.

    """
    switcher = {
        'a0': 0.13,
        'a': 0.21,
        'b': 0.34,
        'c': 0.49,
        'd': 0.76,
    }
    return switcher.get(b_curve, "nothing")


def chi_flex(
             length,
             area,
             MOI,
             f_yield,
             b_curve,
             kapa_BC = None
             ):
    
    # Docstring
    """
    Flexural buckling reduction factor.
    
    Claculates the reduction factor, chi, according to EN1993-1-1 6.3.1.2
    
    Parameters
    ----------
    length : float
        [mm] Column length
    area : float
        [mm^2] Cross section area
    MOI : float
        [mm^4] Moment of inertia
    f_yield : float
        [MPa] Yield stress.
    b_curve : str
        [_] Name of the buckling curve as obtained from Table 6.2 of [1].
        Valid options are {'a0', 'a', 'b', 'c', 'd'}
    kapa_BC : float, optional
        [_] length correction for the effect of the boundary conditions.
        Default = 1, which implies simply supported column
    
    Returns
    -------
    float
        [_] Reduction factor.
    
    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings. Brussels: CEN, 2005.

    """
    if kapa_BC is None:
        kapa_BC = 1.
    
    lmda = lmbda(
        length = length,
        area = area,
        MOI = MOI,
        kapa_BC = kapa_BC,
        E_modulus = None,
        f_yield = f_yield
        )
    
    alpha = imp_factor(b_curve)
    
    phi = (1 + alpha * (lmda - 0.2) + lmda**2) / 2.
    
    chi = 1 / (phi + sqrt(phi**2 - lmda**2))
    
    if chi > 1.:
        chi = 1.
    
    return chi


def N_b_Rd(
           length, 
           area, 
           MOI, 
           f_yield, 
           b_curve, 
           kapa_BC = None, 
           gamma_M1 = None
           ):
    
    # Docstring
    """
    Flexural buckling resistance.
    
    Verifies the resistance of a column against flexural buckling
    according to EN1993-1-1 6.3.1.1.
    
    Parameters
    ----------
    length : float
        [mm] Column length
    area : float
        [mm^2] Cross section area
    MOI : float
        [mm^4] Moment of inertia
    f_yield : float
        [MPa] Yield stress.
    b_curve : str
        [_] Name of the buckling curve as obtained from Table 6.2 of [1].
        Valid options are {'a0', 'a', 'b', 'c', 'd'}
    kapa_BC : float, optional
        [_] Length correction for the effect of the boundary conditions.
        Default = 1, which implies simply supported column
    gamma_M1 : float, optional
        [_] Partial safety factor.
        Default = 1.
    
    Returns
    -------
    float
        [N] Buckling resistance.
    
    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings. Brussels: CEN, 2005.

    """
    if kapa_BC is None:
        kapa_BC = 1.
    
    if gamma_M1 is None:
        gamma_M1 = 1.
    
    chi = chi_flex(length,
        area,
        MOI,
        f_yield,
        b_curve,
        kapa_BC=kapa_BC)
        
    N_b_Rd = area * f_yield * chi / gamma_M1
    
    return N_b_Rd


###### CIRCULAR PLATE  #####
# The following formulas are taken from table 18-3 of:
# W. D. Pilkey, Formulas for stress, strain, and structural matrices. New York: Wiley, 1994.
#
## Hinged perimeter
### Concentrated force applied on circle
#############(not finished, there is an error, to be revised)###########
#
#def circular_plate(a_L, a_1, W_load, r):
#    poisson = 0.3
#    
#    alfa = r / a_L
#    beta = a_1 / a_L
#    
#    C_1 = (3 + poisson) / (1 + poisson) * (1 + beta ** 2) + 2 * beta ** 2 * log(beta)
#    C_2 = (1 - poisson) / (1 + poisson) * (1 - beta ** 2) + 2 * log(beta)
#    C_3 = (3 + poisson) / (1 + poisson) - beta ** 2 * (1 - poisson) / (1 + poisson)
#    
#    if alfa <= beta:
#        M_r = (1 / 4) * W_load * a_L * (1 + poisson) * C_2
#        M_phi = M_r
#    else:
#        M_r = (1 / 4) * (W_load * a_1) / (alfa ** 2) * ((1 - poisson) * (1 - alfa ** 2) * beta ** 2 - 2 * (1 + poisson) * alfa ** 2 * log(alfa))
#        M_phi = (1 / 4) * (W_load * a_1) / (alfa ** 2) * (2 * (1 - poisson) * alfa ** 2 - (1 - poisson) * (1 + alfa ** 2) * beta ** 2 - 2 * (1 + poisson) * alfa ** 2 * log(alfa))
#    
#    return M_r, M_phi