import numpy as np
from itertools import product
from shutil import rmtree
from os import path

# Import the method to be run parametrically
from closed_polygons import cs_calculator as parametric_function


session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

# Give a project name
prj_name = 'closed_polygon'
try:
    prj_name = float(sys.argv[-1])
except:
    pass

# Define a list of lists of input values.
# The parrent list has length of the number of parameters that vary
# The children lists are ranges of values for each parameter.
# A list of all possible combinations is stored in "combinations"
combinations = list(
    product(
        [6], 
	    [30] 
	    )
    )

# Open a file to collect the results
out_file = open('./'+prj_name+'_info.dat', 'a')

# Loop through the combinations of the given input values
for parameter in combinations:
    job_ID = (''.join("%03d-"%e for e in parameter) + prj_name)
    job_ID = job_ID.translate(None, '.')
    if (path.isdir("./" + job_ID)):
        print("Job already exists: A directory with the same name"+job_ID+" exists in the cwd")
        continue
    
    # Make a new subdirectory for the current session
    os.mkdir(job_ID)
    
    # Change working directory
    os.chdir('./'+job_ID)
    
    # The function to be run for the full factorial parametric is called here
    print('Running job: ' + job_ID)
    try:
        job_return = parametric_function(
        n_sides = parameter[0], 
        p_classification = parameter[1], 
        )
    except:
        print('Problem while executing job: '+ job_ID)
        print('Job is canceled and the folder is deleted. See log file (to be writen)')
        Mdb()
        os.chdir('../')
        rmtree(job_ID)
        
    # Return to parent directory
    os.chdir('../')
    
    # Write each returned string to the file separated by newlines
    job_return = str(job_return)
    out_file.write(job_ID + " : " + job_return + "\n")

out_file.close()
