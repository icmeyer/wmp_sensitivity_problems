import numpy as np

import sens_helpers

nmesh = 5
adj_nbatches = 1000
adj_nparticles = 1000
clutch_nbatches = 1000
clutch_nparticles = 1000
# Sensitivities will be evaluated w.r.t. multipole parameters as well
# as across the energy bins defined by sens_e_grid by default
rxns = ['fission', 'absorption', 'elastic']
nuclides = ['Al27', 'U235', 'U238']
sens_e_grid = np.geomspace(1e-5, 2e7, 20)
# Clean xml directory, include '/' at end of directory string
xml_dir = 'heu-sol-therm-020/'

sens_helpers.sequence.run_sensitivity_sequence(xml_dir, nmesh, adj_nbatches,
                                               adj_nparticles, clutch_nbatches,
                                               clutch_nparticles, nuclides, 
                                               rxns, sens_e_grid,
                                               folder_suffix='_{:d}'.format(nmesh))
