from math import sqrt, pi

###### Simply supported plate
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


###### Cylindrical shells
def N_cr_shell(thickness, radius, length):
    # Convert inputs to floats
    thickness, radius, length = float(thickness), float(radius), float(length)
    
    # Elastic critical load acc to EN3-1-6 Annex D
    N_cr_shell = 2 * pi * radius * thickness * sigma_x_Rcr(thickness, radius, length)
    
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
    return sigma_cr
