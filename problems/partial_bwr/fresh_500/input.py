import numpy as np

import sens_helpers

nmesh = 21
adj_nbatches = 100
adj_nparticles = int(1e5)
clutch_nbatches = 100
clutch_nparticles = int(1e5)
# Sensitivities will be evaluated w.r.t. multipole parameters as well
# as across the energy bins defined by sens_e_grid by default
rxns = ['fission', 'absorption', 'elastic']
nuclides = ['H1', 'H2',
            'O16', 'O17', 
            'Fe54', 'Fe56', 'Fe57', 'Fe58',
            'Cr50', 'Cr52', 'Cr53', 'Cr54',
            'Ni58', 'Ni60', 'Ni61', 'Ni62', 'Ni64', 
            'Zr90', 'Zr91', 'Zr92', 'Zr94', 'Zr96',
            'Gd152','Gd154', 'Gd155', 'Gd156', 'Gd157', 'Gd158', 'Gd160',
            'U234', 'U235', 'U236', 'U238']
sens_e_grid = np.geomspace(1e-5, 2e7, 500)
# Clean xml directory, include '/' at end of directory string
xml_dir = 'clean/'

sens_helpers.sequence.run_sensitivity_sequence(xml_dir, nmesh, adj_nbatches,
                                               adj_nparticles, clutch_nbatches,
                                               clutch_nparticles, nuclides, 
                                               rxns, sens_e_grid)
