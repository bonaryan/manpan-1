# Python method converting a 3D cell array from a matlab file (.mat) to an equivalent pickled object containing 3D nested lists.
# scipy is not available in the Python version Abaqus is using
# This script has to be executed by an generic Python of the same version as the one in Abaqus (i.e 2.7.3). Requires pickle and scipy
# The filename is given without the extension

# Module imports
import scipy.io as sio
import pickle

def mat2pkl(filename):

# Load the matlab file
    database = sio.loadmat(filename+'.mat')

# Scipy imports the data of the .mat in a dictionary.
# Get the lists from inside the dictionary
    ppp = database[filename]

# Export with pickle to a .pkl file
    pickle.dump(ppp, open( filename+".pkl", "wb" ))