import numpy as np
from itertools import product
from shutil import rmtree
from os import path

#### INPUT ####
# Import the method to be run parametrically
import closed_polygons

# Give a project name
prj_name = 'imp-I'

# Remove folders after execution? (keep only the output file)
remove_folders = False

# Define a list of lists of input values.
# The parrent list has length of the number of parameters that vary
# The children lists are ranges of values for each parameter.
# A list of all possible combinations is stored in "combinations"
combinations = list(
    product(
        [12], 
        [42],
        ['fcA'],
        [2, 4]
        )
    )

#### END INPUT ####

# Open a file to collect the results
out_file = open('./'+prj_name+'_info.dat', 'a')

# Loop through the combinations of the given input values
for parameter in combinations:
    #job_ID = (''.join("%03d-"%e for e in parameter) + prj_name)
    job_ID = ("%03d-%03d-"%(parameter[0], parameter[1]) + parameter[2] + '-' + "%03d-"%parameter[3] + prj_name)
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
        job_return = closed_polygons.modeler(
            n_sides = parameter[0], 
            p_classification = parameter[1], 
            n_of_waves = parameter[0] / parameter[3], 
            m_of_waves = parameter[0] / parameter[3], 
            fab_class = parameter[2], 
            f_yield = 760, 
            nominal_fy = 'S650',
            IDstring = job_ID, 
            proj_name = prj_name, 
            )
        
        # Return to parent directory
        os.chdir('../')
        
        # Write each returned string to the file separated by newlines
        job_return = str(job_return)
        out_file.write(job_ID + ", " + job_return + "\n")
        
    except:
        print('Problem while executing job: '+ job_ID)
        print('Job is canceled. See log file (no log file yet)')
        os.chdir('../')
    
    # Remove job's folder (only the output information is kept)
    if remove_folders is True:
        rmtree(job_ID)

out_file.close()
