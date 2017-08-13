import numpy as np
from itertools import product
from os import path
from closed_polygons import single_closed_polygon as model_builder


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
        range(5, 7), 
	    [24] 
	    #[1, 2, 3], 
	    #[1, 2, 3]
	    )
    )



# Loop through the combinations of the given input values
for model_params in combinations:
    imp_modes = list(divisorGenerator(model_params[0]))[1:-1]
    for mode in imp_modes:
        spillover = model_params[0] / mode
        IDstring = (''.join("%03d-"%e for e in model_params) + "%03d-"%spillover + prj_name)
        if (path.isdir("./" + IDstring)):
            print("Job already exists: A directory with the same name"+IDstring+" exists in the cwd")
            continue
        
        # A modelling script is called here.
        # The number of given inputs must correspond to the 
        # number of lists in the given combinations
        print('Running model ' + IDstring)
        model_builder(model_params[0], model_params[1], spillover, IDstring)
