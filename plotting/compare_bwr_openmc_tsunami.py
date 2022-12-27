import openmc
import os

import sens_helpers
import scalepy

def hist_y(y):
    return np.hstack([y[0], y])

def plot_sens_err(ax, e_bins, sens, sens_sd, color='b', label=''):
    ax.step(e_bins, hist_y(sens), sens_sd, label=label, color=color)
    midpoints = (e_bins[:-1] + e_bins[1:])/2
    ax.errorbar(midpoints, sens, sens_sd, color=color, fmt=' ', capsize=1)


figs_folder = 'openmc_tsunami_compare/'
os.system('mkdir -p {:s}'.format(figs_folder))

# Import sensitivities from tsunami sdf
sdf_path = '/home/icmeyer/research/gd_pincell/tsunami3d/iso_rings/gd_pincells.sdf'
tsunami_sens_dict = scalepy.tsunami.import_sdf(sdf_path)

# Import sensitivities from openmc
sens_dir = '/home/icmeyer/research/wmp_sensitivity_problems/problems/partial_bwr/fresh/sensitivity_run/'
tally_xml = sens_dir + 'tallies.xml'
tally_out = sens_dir + 'tallies.out'

openmc_sens_dict = sens_helpers.plotting.get_all_tallies(tally_xml, tally_out)

nucs = list(set(openmc_sens_dict.keys()))

for nuc in nucs:
    if nuc in tsunami_sens_dict:
        for reaction in openmc_sens_dict[nuc]['cross_section']:
            print(reaction)
            input()
            e_bins = openmc_sens_dict[nuc][var][reaction]['energy']
            sens = openmc_sens_dict[nuc][var][reaction]['mean']
            sd = openmc_sens_dict[nuc][var][reaction]['sd']


