# Homology modeling with multiple templates
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class
log.verbose()    # request verbose outputenv = environ()  # create a new MODELLER environment to build this model in
# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']
a = automodel(env,
              alnfile  = 'align-multiple.ali', # alignment
 filename
              knowns   = ('1CDO'),     # codes of
 the templates
              sequence = 'adh1DANRER')               # code of the
 target
a.starting_model= 1                 # index of the first model
a.ending_model  = 5                 # index of the last model
# (determines how many models to calculate)
a.make()                            # do the actual homology modeling
