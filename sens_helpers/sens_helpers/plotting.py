import pickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import openmc
import csv
import xml.etree.ElementTree as ET
import os
import collections
import pickle
import re
import shutil

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

def plot_slice(fig, ax, data, axis, n, nmesh, error=False):
    if np.any(~np.isnan(data)):
        vmax  = np.max(data[~np.isnan(data)])
        vmin  = np.min(data[~np.isnan(data)])
    else:
        vmax = 0
        vmin = 0

    cmap = 'jet'
    gcf = ax.imshow(data, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(gcf, cax=cax)
    if error:
        cbar.set_label('% Error')
    else:
        cbar.set_label('')
    
    ax.set_title('{:s}-Axis Slice: {:d}/{:d}'.format(axis, n+1, nmesh))

def plot_source():
    with open('source.pkl','rb') as f:
        source_dict = pickle.load(f)
    source = source_dict['source']
    source_sd = source_dict['source_sd']

    nmesh = int(np.round(len(source)**(1/3)))
    source = source.reshape([nmesh, nmesh, nmesh], order='C')
    source_sd = source_sd.reshape([nmesh, nmesh, nmesh], order='C')

    figures = []
    for n in range(nmesh):
        fig = plt.figure(figsize=(11,8))
        fig.suptitle('Source Tally')

        means = [source[n, :, :], source[:, n, :], source[:, :, n]]
        sds = [source_sd[n, :, :], source_sd[:, n, :], source_sd[:, :, n]]
        axis_names = ['z', 'y', 'x']

        for i in range(3):
            ax = fig.add_subplot(2,3,i+1)
            plot_slice(fig, ax, means[i], axis_names[i], n, nmesh)
            ax = fig.add_subplot(2,3,i+4)
            plot_slice(fig, ax, means[i]/sds[i]*100, axis_names[i], n, nmesh, error=True)

        fig.tight_layout()
        figures.append(fig)

    return figures

def plot_adjoint():
    with open('adjoint.pkl','rb') as f:
        adjoint = pickle.load(f)

    nmesh = int(np.round(len(adjoint)**(1/3)))
    adjoint = adjoint.reshape([nmesh, nmesh, nmesh], order='C')

    figures = []
    for n in range(nmesh):
        fig = plt.figure(figsize=(11,4))
        fig.suptitle('Adjoint')

        means = [adjoint[n, :, :], adjoint[:, n, :], adjoint[:, :, n]]
        axis_names = ['z', 'y', 'x']

        for i in range(3):
            ax = fig.add_subplot(1,3,i+1)
            plot_slice(fig, ax, means[i], axis_names[i], n, nmesh)

        fig.tight_layout()
        figures.append(fig)

    return figures

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def get_tally_result(tallies_result, info):
    vals = []
    errs = []
    regex_comp = re.compile("TALLY " + str(info['tallyid']))
    found = False
    with open(tallies_result) as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if found:
                if line[:2]==' =':
                    table_end = i-1
                    break
                if i>=(table_start+5):
                    elements = line.split()
                    vals.append(elements[2])
                    errs.append(elements[4])
            else:
                for match in re.finditer(regex_comp, line):
                    found=True
                    table_start=i
                    continue
    vals = np.array(vals, dtype=float)
    errs = np.array(errs, dtype=float)
    result = {'mean': vals,
              'sd': errs}
    return result

def nested_dict():
    return collections.defaultdict(nested_dict)

def get_tally(sens):
    tallyid = sens.attrib['id']
    nuclide = sens.find('nuclide').text
    variety = sens.find('variable').text

    reaction = sens.find('reaction')
    if reaction is None:
        reaction = 'none'
    else:
        reaction = reaction.text

    energy = sens.find('energy')
    if energy is not None:
        energy = np.fromstring(sens.find('energy').text, dtype=float, sep=' ')

    return {'nuclide': nuclide,
            'var': variety,
            'reaction': reaction,
            'energy': energy,
            'tallyid': tallyid}

def get_all_tallies(tally_xml, tally_out):
    tree = ET.parse(tally_xml)
    root = tree.getroot()
    sens_lib = nested_dict()
    for sens in root.findall('sensitivity'):
        info = get_tally(sens)
        result = get_tally_result(tally_out, info)
        result['energy'] = info['energy']
        sens_lib[info['nuclide']][info['var']][info['reaction']] = result
        print(info['nuclide'], info['var'], info['reaction'])
    return sens_lib

def hist_y(x):
    return np.hstack([x[0], x])

def plot_xs_sens(sens_dict, nuc, var, reaction, title):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    print(nuc, var, reaction)
    e_bins = sens_dict[nuc][var][reaction]['energy']
    sens = sens_dict[nuc][var][reaction]['mean']
    sd = sens_dict[nuc][var][reaction]['sd']
    mids = (e_bins[:-1] + e_bins[1:])/2
    step = ax.step(e_bins, hist_y(sens))
    ax.errorbar(mids, sens, yerr=sd, capsize=2, capthick=1, fmt=' ', color=step[0].get_color())
    ax.set_ylabel('Sensitivity')
    ax.set_xlabel('Energy (eV)')
    ax.set_title(title)
    ax.set_xscale('log')
    ax.set_yscale('symlog')
    fig.tight_layout()
    return fig

def plot_wmp_sens(sens_dict, nuc, var, reaction, title):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    print(nuc, var, reaction)
    sens = sens_dict[nuc][var]['none']['mean']
    x = np.arange(len(sens))
    sd = sens_dict[nuc][var]['none']['sd']
    ax.errorbar(x, sens, yerr=sd, capsize=2, capthick=1, fmt='.')
    ax.set_ylabel('Sensitivity')
    ax.set_xlabel('Energy (eV)')
    ax.set_title(title)
    ax.set_yscale('symlog')
    fig.tight_layout()
    return fig

def plot_all_sens(nuclides, reactions, var_words):
    figs = []
    tally_xml = 'tallies.xml'
    tally_out = 'tallies.out'
    sens_dict = get_all_tallies(tally_xml, tally_out)
    for nuc in nuclides:
        for reaction in reactions:
            for var_word in var_words:
                title = nuc+' '+reaction+' '+var_word
                if var_word in ['multipole', 'curve fit']:
                    fig = plot_wmp_sens(sens_dict, nuc, var_word, reaction, title)
                elif var_word == 'cross_section':
                    fig = plot_xs_sens(sens_dict, nuc, var_word, reaction, title)
                figs.append(fig)
    return figs


def plot_sens_test():
    rxns = ['fission', 'absorption', 'elastic']
    var_words = ['cross_section', 'multipole', 'curve_fit']
    nuclides = ['Al27', 'U235', 'U238']
    plot_all_sens(nuclides, rxns, var_words)


if __name__=='__main__':
    import os
    sens_dir = 'heu-sol-therm-020_nmesh_11'
    os.chdir(sens_dir)
    # plot_source()
    # plot_adjoint()
    plot_sens_test()
    plt.show()

    