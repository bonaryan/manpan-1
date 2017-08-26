import numpy as np
from itertools import product
from shutil import rmtree
from os import path

#### INPUT ####
# Import the method to be run parametrically
import closed_polygons
import abq_toolset
import EN_tools

# Give a project name
prj_name = 'specimens'

# Remove folders after execution? (keep only the output file)
remove_folders = True

# Define a list of lists of input values.
# The parrent list has length of the number of parameters that vary
# The children lists are ranges of values for each parameter.
# A list of all possible combinations is stored in "combinations"
combinations = list(
    product(
        [16, 20, 24], 
        [42]
        )
    )

#### END INPUT ####

# Open a file to collect the results
out_file = open('./'+prj_name+'_info.dat', 'a')

# Loop through the combinations of the given input values
for parameter in combinations:
    job_ID = (''.join("%03d-"%e for e in parameter) + prj_name)
    #job_ID = ("%03d-%03d-"%(parameter[0], parameter[1]) + parameter[2] + '-' + "%03d-"%parameter[3] + prj_name)
    job_ID = job_ID.translate(None, '.')
    if (path.isdir("./" + job_ID)):
        print("Job already exists: A directory with the same name"+job_ID+" exists in the cwd")
        continue
    
    # Make a new subdirectory for the current session
    os.mkdir(job_ID)
    
    # Change working directory
    os.chdir('./' + job_ID)
    
    # The function to be run for the full factorial parametric is called here
    print('Running job: ' + job_ID)
    
    try:
        thickness = 2.
        r_circle = closed_polygons.class_2_radius(
            n_sides = parameter[0],
            p_classification = parameter[1],
            thickness = thickness,
            f_yield = 650.
            )
        
        profile = closed_polygons.cs_calculator(
            n_sides = parameter[0], 
            r_circle = r_circle, 
            p_classification = parameter[1], 
            column_length = 1000.,
            f_yield = 650., 
            )
        
        x_corners = profile[-2]
        y_corners = profile[-1]
        
        nodes = [x_corners, y_corners]
        elem = [range(0, len(x_corners)), range(1, len(x_corners)) + [0], len(x_corners) * [thickness]]
        cs_properties = abq_toolset.cs_prop(nodes, elem)
        area = cs_properties[0]
        I_2 = cs_properties[-2]
        
        lmda = EN_tools.lmbda(
            1000.,
            area,
            I_2,
            kapa_BC = 1.,
            E_modulus = 210000.,
            f_yield = 650.
            )
        
        job_return = [r_circle]+[list(profile[0:3])+[lmda]]
    
    except:
        print('Problem while executing job: '+ job_ID)
        print('Job is canceled. See log file (no log file yet)')
        os.chdir('../')
    
    # Remove job's folder (only the output information is kept)
    #if remove_folders is True:
    #    rmtree(job_ID)

out_file.close()
