import numpy as np

import sens_helpers

nmesh = 5
adj_nbatches = 100
adj_nparticles = 100
clutch_nbatches = 100
clutch_nparticles = 100
# Sensitivities will be evaluated w.r.t. multipole parameters as well
# as across the energy bins defined by sens_e_grid by default
rxns = ['fission', 'absorption', 'elastic']
nuclides = ['U235', 
            'Gd154',
            'Gd157']
sens_e_grid = np.geomspace(1e-5, 2e7, 20)
# Clean xml directory, include '/' at end of directory string
xml_dir = 'clean/'

sens_helpers.sequence.run_sensitivity_sequence(xml_dir, nmesh, adj_nbatches,
                                               adj_nparticles, clutch_nbatches,
                                               clutch_nparticles, nuclides, 
                                               rxns, sens_e_grid)
