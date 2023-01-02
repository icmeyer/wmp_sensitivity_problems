import openmc
import numpy as np
import scipy
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os
import pickle

import mpullse
import scalepy
import sens_helpers

# sens_folder = '/home/icmeyer/research/wmp_sensitivity_problems/problems/pu-sol-therm-034-015/sensitivity_run'

# OpenMC Sensitivity
sens_dir = '/home/icmeyer/research/wmp_sensitivity_problems/problems/partial_bwr/fresh/sensitivity_run/'

bench_name = sens_dir.split('/')[-3]

figs_dir = bench_name + '_figs/'
os.system('mkdir -p {:s}'.format(figs_dir))

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
    # wmp_obj = mpullse.wmp_data.get_wmp_elemA(nuc, mask_str='coefs')
    wmp_obj = mpullse.wmp_data.get_wmp_elemA(nuc)
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
    dmg_dp = wmp_obj.multigroup_derivative(mp_ebins, flux='overE', 
                                           indexing='rgxp')
    mg = wmp_obj.multigroup(mp_ebins, flux='overE', indexing='rg')

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
    reaction_labels = [r'$\sigma_{s}$', r'$\sigma_{\gamma}$', r'$\sigma_{f}$']
    k_sens_from_mg = []
    k_sens_from_mg_sd = []
    for reaction in reactions:
        sens = openmc_sens_dict[nuc]['cross_section'][reaction]['mean'][el:eh-1]
        sd = openmc_sens_dict[nuc]['cross_section'][reaction]['sd'][el:eh-1]
        k_sens_from_mg.append(sens)
        k_sens_from_mg_sd.append(sd)

    dk_ds_vec = (np.hstack(k_sens_from_mg) / mg).reshape([1, -1])

    # # "Cheating" to get dmg_dp from known sensitivities
    # # proves the math is correct
    # dmg_dp = np.linalg.pinv(dk_ds_vec.T @ dk_ds_vec) @ dk_ds_vec.T @ dk_dp.reshape([1, -1])

    # Calculate eigenvalue derivative from wmp sens
    A = scipy.linalg.pinv(dmg_dp @ dmg_dp.T)
    print(dk_dp.shape, dmg_dp.T.shape, A.shape)
    dk_dmg_from_wmp = dk_dp @ dmg_dp.T @ A
    k_sens_from_wmp = dk_dmg_from_wmp * mg
    k_sens_from_wmp = np.split(k_sens_from_wmp, nrxns)


    dk_dp_from_xs_sens = dk_ds_vec @ dmg_dp

    # Plotting 

    # # Plot multipole sens values
    # fig = plt.figure()
    # ax = fig.add_subplot(311)
    # print('wmp sens shape', wmp_sens.shape)
    # # ax.scatter(np.arange(len(wmp_sens)), wmp_sens)
    # ax.scatter(np.arange(len(wmp_sens)), dk_dp, alpha=0.5)
    # ax.scatter(np.arange(len(wmp_sens)), dk_dp_from_xs_sens, alpha=0.5)
    # ax.set_yscale('symlog')
    # ax = fig.add_subplot(312, sharex=ax)
    # ax.scatter(np.arange(len(wmp_sens)), 
    #                      (dk_dp_from_xs_sens-dk_dp)/dk_dp*100)
    # ax.set_ylabel('Relative Difference [\%]')
    # ax.set_yscale('symlog')
    # ax = fig.add_subplot(313, sharex=ax)
    # ax.scatter(np.arange(len(wmp_sens)), mp_mean)
    # ax.set_yscale('symlog')
    

    fig = plt.figure(figsize=(14,8))
    if fissionable:
        axd = fig.subplot_mosaic(
            """
            aaddg
            bbeeh
            ccffi
            """
        )
        axes = [axd['a'], axd['b'], axd['c']]
        axes2 = [axd['d'], axd['e'], axd['f']]
        leg_axes = [axd['g'], axd['h'], axd['i']]
    else:
        axd = fig.subplot_mosaic(
            """
            aaddg
            bbeeh
            """
        )
        axes = [axd['a'], axd['b']]
        axes2 = [axd['d'], axd['e']]
        leg_axes = [axd['g'], axd['h']]

    edges = wmp_obj.get_window_edges()
    alphas = 1/(np.arange(len(edges))/2 + 1)
    for r, reaction in enumerate(reactions):
        ax = axes[r]
        label = 'OpenMC Sensitivity Tally'
        sens_helpers.plotting.plot_sens_err(ax, mp_ebins,
                                            k_sens_from_mg[r],
                                            k_sens_from_mg_sd[r],
                                            label=label)
        label = 'Derived from WMP'
        sens_helpers.plotting.plot_sens_err(ax, mp_ebins,
                                            k_sens_from_wmp[r],
                                            label=label)

        ax.axvline(x=edges[0], alpha=alphas[0], linestyle='--', color='k', label='WMP Window Edges')
        for i in range(len(edges)):
            ax.axvline(x=edges[i], alpha=alphas[i], linestyle='--', color='k')

        ax.set_xscale('log')
        ax.set_yscale('symlog')
        ax.axhline(y=0, color='k', linewidth=0.5)
        #ax.set_ylabel(r'Sensitivity ($\frac{\sigma_g}{k} \frac{\partial k}{\sigma_g})$')
        ax.set_ylabel(r'Sensitivity')
        if r==(len(reactions)-1):
            ax.set_xlabel(r'Energy [eV]')
        ax.set_title(reaction_labels[r])
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        h, l = ax.get_legend_handles_labels()

        ax = axes2[r]
        sens_helpers.plotting.plot_sens_err(ax, mp_ebins,
                                            k_sens_from_mg[r]/log_width,
                                            (k_sens_from_mg_sd[r]**2/log_width)**(1/2))
        sens_helpers.plotting.plot_sens_err(ax, mp_ebins,
                                            k_sens_from_wmp[r]/log_width)
        for i in range(len(edges)):
            ax.axvline(x=edges[i], alpha=alphas[i], linestyle='--', color='k')
        ax.set_xscale('log')
        ax.set_yscale('symlog')
        ax.set_ylabel('Sensitivity \nper unit lethargy')
        if r==(len(reactions)-1):
            ax.set_xlabel(r'Energy [eV]')
        ax.set_title(reaction_labels[r])
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax.axhline(y=0, color='k', linewidth=0.5)



        ax = leg_axes[r]
        ax.axis('off')
        if r==0:
            ax.legend(h, l, loc='center right')

    fig.suptitle(scalepy.isotopes.name_conversions.ElemA_to_latex(nuc))
    fig.tight_layout()
    fig_title = '{:s}'.format(nuc)

    print(figs_dir + fig_title)
    fig.savefig(figs_dir + fig_title + '.png')
    fig.savefig(figs_dir + fig_title + '.pdf', dpi=600)

    # Plot log energy grid
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # plt.step(mp_ebins[1:], -1/log_width)
    # ax.set_xscale('log')
    
    
    plt.close(fig)



    


