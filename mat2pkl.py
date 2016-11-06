# Python script converting a 3D cell array from a matlab file (.mat) to an equivalent pickled object containing 3D nested lists.
# It is written as an independent script because scipy is required which is not available in the Python version comming with Abaqus
# This script has to be executed by an generic Python with scipy.

# Module imports
import scipy.io as sio
import pickle

# Load the matlab file
database = sio.loadmat('profiles.mat')

# Scipy imports the data of the .mat in a dictionary.
# Get the lists from inside the dictionary
ppp = database['profiles']

# Export with pickle to a .pkl file
pickle.dump(ppp, open( "profiles.pkl", "wb" ))