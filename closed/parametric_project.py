import numpy as np
from itertools import product
from shutil import rmtree
from os import path

# Import the method to be run parametrically
from closed_polygons import single_closed_polygon as parametric_function


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

# Loop through the combinations of the given input values
for parameter in combinations:
    IDstring = (''.join("%03d-"%e for e in parameter) + prj_name)
    IDstring = IDstring.translate(None, '.')
    if (path.isdir("./" + IDstring)):
        print("Job already exists: A directory with the same name"+IDstring+" exists in the cwd")
        continue
    
    # Make a new subdirectory for the current session
    os.mkdir(IDstring)
    
    # Change working directory
    os.chdir('./'+IDstring)
    
    # The function to be run for the full factorial parametric is called here
    print('Running job: ' + IDstring)
    try:
        parametric_function(parameter[0], parameter[1], parameter[0]/2, IDstring, prj_name)
    except:
        print('Problem while executing job: '+ IDstring)
        print('Job is canceled and the folder is deleted. See log file')
        os.chdir('../')
        rmtree(IDstring)
        
    
    # Return to parent directory
    os.chdir('../')