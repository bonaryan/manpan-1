import numpy as np
from itertools import product
from shutil import rmtree
from os import path
from time import sleep
from random import random

#### INPUT ####
# Import the method to be run parametrically
import EN_tools

# Give a project name
prj_name = 'shell_critical_load'

# Remove folders after execution? (keep only the output file)
remove_folders = True

# Delay jobs to avoid collisions of parallel jobs when running on the cluster
delay_jobs = False

# Define a list of lists of input values.
# The parrent list has length of the number of parameters that vary
# The children lists are ranges of values for each parameter.
# A list of all possible combinations is stored in "combinations"
combinations = list(
    product(
        range(6, 31, 2),
        range(27, 55, 3),
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
    
    # Wait some seconds to avoid multiple initiation of jobs
    if delay_jobs is True:
        sleep(round(10*random(), 2))
    
    # Check if the directory exists
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
        diameter = 500.
        radius = diameter / 2
        epsilon = sqrt(235./ 381.)
        thickness = pi * diameter / (parameter[0] * parameter[1] * epsilon)
        width = pi * diameter / parameter[0]
        #job_return = pi * diameter * thickness * EN_tools.sigma_cr_plate(thickness, width)
        
        # Shell critical load
        job_return = EN_tools.sigma_x_Rcr(thickness, radius, 2*pi*radius)
        N_cr_shell = job_return[0] * thickness * pi * diameter
        
        # Return to parent directory
        os.chdir('../')
        
        # Write each returned string to the file separated by newlines
        job_return = str(job_return) + str(N_cr_shell)
        out_file.write(job_ID + ", " + job_return + "\n")
        
    except:
        print('Problem while executing job: '+ job_ID)
        print('Job is canceled. See log file (no log file yet)')
        os.chdir('../')
    
    # Remove job's folder (only the output information is kept)
    if remove_folders is True:
        rmtree(job_ID)

out_file.close()
