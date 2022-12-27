import openmc
import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
import pickle

import mpullse
import sens_helpers

# sens_folder = '/home/icmeyer/research/wmp_sensitivity_problems/problems/pu-sol-therm-034-015/sensitivity_run'

# OpenMC Sensitivity
sens_dir = '/home/icmeyer/research/wmp_sensitivity_problems/problems/partial_bwr/fresh/sensitivity_run/'

bench_name = sens_dir.split('/')[-2]

# pickle handling
pickle_folder = 'pickles/'
os.system('mkdir -p {:s}'.format(pickle_folder))
omc_pkl_path = pickle_folder + 'omc_sens_{:s}.pkl'.format(bench_name)
if os.path.exists(omc_pkl_path):
    print('Loading sensitivities from pickle')
    with open(omc_pkl_path, 'rb') as f:
        openmc_sens_dict = pickle.load(f)
else:
    tally_xml = sens_dir + 'tallies.xml'
    tally_out = sens_dir + 'tallies.out'
    openmc_sens_dict = sens_helpers.plotting.get_all_tallies(tally_xml, tally_out)
    with open(omc_pkl_path, 'wb') as f:
        pickle.dump(openmc_sens_dict, f)

# Get list of nuclides
nucs = np.sort(list(set(openmc_sens_dict.keys())))
# nucs = ['H1', 'Gd155']
# nucs = ['H1']

# Get energy structure of xs sens
e_bins = openmc_sens_dict[nucs[0]]['cross_section']['absorption']['energy']

for nuc in nucs:
    print(nuc)
    # wmp_obj = mpullse.wmp_data.get_wmp_elemA(nuc, mask_str='poles')
    wmp_obj = mpullse.wmp_data.get_wmp_elemA(nuc, mask_str='coefs')
    # wmp_obj = mpullse.wmp_data.get_wmp_elemA(nuc)
    mp_mean = wmp_obj.get_par_vector()
    fissionable = wmp_obj.fissionable
    nrxns = 2 + fissionable

    # Find indices of wmp energy range
    el = np.searchsorted(e_bins, wmp_obj.E_min)
    eh = np.searchsorted(e_bins, wmp_obj.E_max)
    print(wmp_obj.E_min, wmp_obj.E_max)

    mp_ebins = e_bins[el:eh]
    log_width = np.log(mp_ebins[1:]) - np.log(mp_ebins[:-1])

    # Calculate xs deriv wrt wmp (dsigma/dpi)
    dmg_dp = wmp_obj.multigroup_derivative(mp_ebins, flux='flat', 
                                           indexing='rgxp')
    mg = wmp_obj.multigroup(mp_ebins, flux='flat', indexing='rg')

    # Get wmp sensitivity from OpenMC sensitivity run
    mp_sens = openmc_sens_dict[nuc]['multipole']['none']['mean']
    cf_sens = openmc_sens_dict[nuc]['curve_fit']['none']['mean']
    print(openmc_sens_dict[nuc].keys())
    print('poles', mp_sens.shape, 'curve_fit', cf_sens.shape)
    wmp_sens = np.hstack([mp_sens, cf_sens])
    if wmp_obj.is_masked:
        wmp_sens = wmp_sens[wmp_obj.cov_mask]

    print('N parameters:', wmp_obj.nparams)
    print('N mp sensitivities:', wmp_sens.shape)
    dk_dp = wmp_sens # * mp_mean
    # dk_dp = wmp_sens 


    # Get mg sens from sensitivity run
    reactions = ['elastic', 'absorption'] + fissionable*['fission']
    k_sens_from_mg = []
    k_sens_from_mg_sd = []
    for reaction in reactions:
        sens = openmc_sens_dict[nuc]['cross_section'][reaction]['mean'][el:eh-1]
        sd = openmc_sens_dict[nuc]['cross_section'][reaction]['sd'][el:eh-1]
        k_sens_from_mg.append(sens)
        k_sens_from_mg_sd.append(sd)

    dk_ds_vec = (np.hstack(k_sens_from_mg) / mg).reshape([1, -1])

    # # "Cheating" to get dmg_dp from known sensitivities
    # dmg_dp = np.linalg.pinv(dk_ds_vec.T @ dk_ds_vec) @ dk_ds_vec.T @ dk_dp.reshape([1, -1])

    # Calculate eigenvalue derivative from wmp sens
    A = scipy.linalg.pinv(dmg_dp @ dmg_dp.T)
    print(dk_dp.shape, dmg_dp.T.shape, A.shape)
    dk_dmg_from_wmp = dk_dp @ dmg_dp.T @ A
    k_sens_from_wmp = dk_dmg_from_wmp * mg
    k_sens_from_wmp = np.split(k_sens_from_wmp, nrxns)


    dk_dp_from_xs_sens = dk_ds_vec @ dmg_dp


    # Plotting 
    fig = plt.figure()
    ax = fig.add_subplot(311)
    print('wmp sens shape', wmp_sens.shape)
    # ax.scatter(np.arange(len(wmp_sens)), wmp_sens)
    ax.scatter(np.arange(len(wmp_sens)), dk_dp, alpha=0.5)
    ax.scatter(np.arange(len(wmp_sens)), dk_dp_from_xs_sens, alpha=0.5)
    ax.set_yscale('symlog')
    ax = fig.add_subplot(312, sharex=ax)
    ax.scatter(np.arange(len(wmp_sens)), 
                         (dk_dp_from_xs_sens-dk_dp)/dk_dp*100)
    ax.set_ylabel('Relative Difference [\%]')
    ax.set_yscale('symlog')
    ax = fig.add_subplot(313, sharex=ax)
    ax.scatter(np.arange(len(wmp_sens)), mp_mean)
    ax.set_yscale('symlog')
    

    fig = plt.figure()
    if fissionable:
        axd = fig.subplot_mosaic(
            """
            a
            b
            c
            """
        )
        axes = [axd['a'], axd['b'], axd['c']]
    else:
        axd = fig.subplot_mosaic(
            """
            a
            b
            """
        )
        axes = [axd['a'], axd['b']]

    for r, reaction in enumerate(reactions):
        ax = axes[r]
        label = 'Cross Section Tally'
        sens_helpers.plotting.plot_sens_err(ax, mp_ebins,
                                            k_sens_from_mg[r]/log_width,
                                            k_sens_from_mg_sd[r],
                                            label=label)
        label = 'Derived from WMP'
        sens_helpers.plotting.plot_sens_err(ax, mp_ebins,
                                            k_sens_from_wmp[r]/log_width,
                                            label=label)
        ax.set_xscale('log')
        ax.set_yscale('symlog')
        ax.set_title(reaction)
        ax.legend()
    fig.suptitle(nuc)
    plt.show()
    plt.close(fig)



    


