import numpy as np
import abq_toolset as xtr
import os
import sys
import odbAccess

out_file1 = open('./batch_info.dat', 'a')

for n_sides in range(5, 26 + 1):
    for p_classification in [24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54]:
        
        # GEOMETRY AND MATERIAL
        
        # Number of requested eigenvalues
        n_eigen = 20
        
        # Yield stress
        f_yield = 381.
        
        # GEOMETRY
        r_circle = 250
        
        # Radius of the polygons circumscribed
        r_circum = (pi * r_circle) / (n_sides * sin(pi / n_sides))
        
        # Diameter
        diameter = 2 * r_circum
        
        # Column length
        column_length = 2 * diameter
        
        # Central angles
        theta = 2 * pi / n_sides
        
        # Width of each side
        w_side = diameter * sin(pi / n_sides)
        
        # Perimeter
        perimeter = n_sides * diameter * sin(theta / 2)
        
        # Epsilon for the material
        epsilon = sqrt(235. / f_yield)
        
        # Thickness for profile on class 3-4 limit (classification as plated, not tube)
        shell_thickness = (diameter * sin(theta/2)) / (p_classification * epsilon)
        
        # List of angles of the polygon corner points to the x-axis
        phii = []
        for i_index in range(n_sides):
            phii.append(i_index * theta)
        
        # Polygon corners coordinates
        x_corners = r_circum * np.cos(phii)
        y_corners = r_circum * np.sin(phii)
        
        # END GEOMETRY AND MATERIAL
        
        # CS PROPERTIES
        
        # Classification as tube
        t_classification = 2 * r_circle / (shell_thickness * epsilon ** 2)
        
        # Calculate cross-sectional properties using the function from abq_toolset
        coord = [x_corners, y_corners]
        ends = [range(n_sides), range(1,n_sides)+[0], [shell_thickness]*(n_sides)]
        
        Area, xc, yc, Ix, Iy, Ixy, I1, I2, theta_principal = xtr.cs_prop(coord, ends)
        
        # END CS PROPERTIES
        
        # CRITICAL STRESS
        
        # Elastic critical load acc. to EN3-1-5 Annex A
        # Uniform compression is assumed (k_sigma = 4)
        psi = 1.
        kapa_sigma = 8.2 / (1.05 + psi)
        sigma_cr_plate =  190000 * (shell_thickness / w_side) ** 2
        N_cr_plate = Area * kapa_sigma * sigma_cr_plate
        
        # Elastic critical load acc. to EN3-1-6 Annex D
        omega = column_length / sqrt(r_circle * shell_thickness)
        if 1.7 <= omega and omega <= (r_circle / shell_thickness):
            C_x = 1.
            length_category = 'medium'
        elif omega < 1.7:
            C_x = 1.36 - 1.83 / omega + 2.07 / omega ** 2
            length_category = 'short'
        else:
            # C_x_b is read on table D.1 of EN3-1-6 Annex D acc. to BCs
            # BC1 - BC1 is used on the Abaqus models (both ends clamped, see EN3-1-5 table 5.1)
            C_x_b = 6.
            C_x_N = min((1 + 0.2 * (1 - 2 * omega * shell_thickness / r_circle) / C_x_b), 0.6)
            C_x = C_x_N
            length_category = 'long'
        
        # Calculate critical stress, eq. D.2 on EN3-1-5 D.1.2.1-5
        sigma_x_Rcr = 0.605 * 210000 * C_x * shell_thickness / r_circle
        
        # Elastic critical load acc to EN3-1-6 Annex D
        N_cr_shell = Area * sigma_x_Rcr
        
        # END CRITICAL STRESS
        
        # DESIGN RESISTANCE
        
        # Aeff calculation.
        # Reduction factor for the effective area of the profile acc. to EC3-1-5
        lambda_p = p_classification / (28.4 * sqrt(kapa_sigma))
        if lambda_p > 0.673 and int(p_classification) > 42:
            rho = (lambda_p - 0.055 * (3 + psi)) / lambda_p ** 2
        else:
            rho = 1.
        
        # Effective area
        A_eff = rho * Area
        
        # Axial compression resistance , Npl
        N_pl_Rd = A_eff * f_yield
        
        # END DESIGN RESISTANCE
        
        # READ EIGEN VALUES FROM ODB
        # Model filename and directory
        IDstring = "%03d-%03d"%(n_sides, p_classification)
        
        # Enter the subdirectory where the specific model is located
        os.chdir(IDstring)
        
        # Open the buckling odb
        bckl_odb = odbAccess.openOdb(path='BCKL-'+IDstring+'.odb')
        bckl_step = bckl_odb.steps['bckl']
        
        # Gather the eigenvalues from the .odb files in the variable "eigenvalues"
        # and print them in the "eigen_string"
        eigenvalues = ()
        eigen_string = ""
        for J_eigenvalues in range(1, n_eigen + 1):
            current_eigen = float(bckl_step.frames[J_eigenvalues].description[-11:])
            eigenvalues = eigenvalues + (current_eigen, )
            eigen_string = eigen_string + "%.3E "%(current_eigen)
        
        # END READ EIGEN VALUES FROM ODB        
        
        # Return to parent directory
        os.chdir('..')
        
        # Create and populate an output text with model information
        out_file1.write("%03d %07.3f %07.3f %07.3f %07.3f %07.3f %010.3f %03d %05.3f %05.1f %05.1f %05.3f %010.3f %.3E %09.3f %09.3f %.3E %.3E "%(n_sides, 2 * r_circle, 2 * r_circum, w_side, perimeter, shell_thickness, Area, f_yield, epsilon, p_classification, t_classification, rho, A_eff, N_pl_Rd, sigma_cr_plate, sigma_x_Rcr, N_cr_plate, N_cr_shell))
        out_file1.write(eigen_string+"\n")

out_file1.close()
