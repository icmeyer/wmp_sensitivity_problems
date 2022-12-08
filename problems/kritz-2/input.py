import numpy as np

import sens_helpers

nmesh = 10
adj_nbatches = 100
adj_nparticles = 100000
clutch_nbatches = 100
clutch_nparticles = 100000
# Sensitivities will be evaluated w.r.t. multipole parameters as well
# as across the energy bins defined by sens_e_grid by default
rxns = ['fission', 'absorption', 'elastic']
nuclides = 'KEY'
sens_e_grid = sens_helpers.GROUPS['SCALE-252'][::-1]
# Clean xml directory, include '/' at end of directory string
xml_dir = 'clean/'

sens_helpers.sequence.run_sensitivity_sequence(xml_dir, nmesh, adj_nbatches,
                                               adj_nparticles, clutch_nbatches,
                                               clutch_nparticles, nuclides, 
                                               rxns, sens_e_grid)
