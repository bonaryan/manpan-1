# This script should define a function to write on a text file some
# characteristics of the model (examples mentioned below)
# Additionaly, a pickle file can be exported containing all the script variables
# so that when loaded, the objects of the model could be restored and the
# can be modified.

# Info to be written in a text file:
# column lengths, total length 
# plate thicknesses
# area
# N_rd
# number of bolts per span?
# classification as plate and as tube
# resulted chi reduction factor

# Create a text file to write model information
infofile = open("mdl_nifo.txt", "w")
infofile.write('a = '+str(a)\n)

