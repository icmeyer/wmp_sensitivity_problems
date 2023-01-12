import numpy as np
import shutil
import h5py
import copy
import os
import pickle
import re
import socket
from pathlib import Path

import openmc

# Getting template files
from .data import KEY_NUCLIDES
from .openmc_help import get_nucs_from_mat_xml
import importlib.resources as pkg_resources
from . import templates  
tally_template = pkg_resources.read_text(templates, 'tallies.xml')



def mkdir(directory):
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.mkdir(directory)

def replace_line(filename, old, new_str):
    new = ""
    regex = re.compile(old)
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            line = re.sub(regex, new_str, line)
            print(line)
            new += line
            new += "\n"
    with open(filename, 'w') as f:
        f.write(new)

def initialize_adjoint_run(clean_dir, nmesh, nbatches, nparticles):
    # import the current simulation and output an appropriate tally file 
    # to initiate the simulation
    print('loading file', clean_dir + 'materials.xml')
    materials = openmc.Materials.from_xml(clean_dir + 'materials.xml')
    settings = openmc.Settings.from_xml(clean_dir + 'settings.xml')
    geometry = openmc.Geometry.from_xml(clean_dir + 'geometry.xml')

    settings.temperature['method'] = 'nearest'
    settings.temperature['tolerance'] = 1000
    settings.temperature['multipole'] = True

    settings.particles = nparticles
    settings.batches = nbatches
    settings.inactive = 20

    # Make mesh for problem
    box = geometry.bounding_box
    lengths = box[1] - box[0]
    pitch = lengths/nmesh
    lattice = openmc.RectLattice()

    # Determine if problem is 2D
    if np.isinf(box[1][2]):
        is2D = True
    else:
        is2D = False

    if is2D:
        lattice.lower_left = box[0][:2]
        lattice.pitch = pitch[:2]
    else:
        lattice.lower_left = box[0]
        lattice.pitch = pitch


    # Fill lattice with universes
    u = openmc.Universe()
    if is2D:
        univ_array = np.empty((nmesh, nmesh), dtype=openmc.Universe)
        univ_array[:, :] = u
    else:
        univ_array = np.empty((nmesh, nmesh, nmesh), dtype=openmc.Universe)
        univ_array[:, :, :] = u
    lattice.universes = univ_array
    # Extract mesh
    mesh =  openmc.RegularMesh.from_rect_lattice(lattice)


    # Export all the xmls here
    print('Exporting adjoint run xmls')
    materials.export_to_xml()
    settings.export_to_xml()
    geometry.export_to_xml()

    # tallies file currently needs to be manually edited xml
    # remove current if present and copy clean version
    cur_tally = './tallies.xml'
    if os.path.exists(cur_tally):
        os.remove(cur_tally)
    with open(cur_tally, 'w') as f:
        for line in tally_template:
            f.write(line)

    # Write mesh info to tally xml
    old = "MESH_XML"
    mesh_lines = ''
    if is2D:
        mesh_lines += '  <dimension>{:d} {:d} </dimension>'\
                      '\n'.format(nmesh,nmesh)
        mesh_lines += '  <lower_left>{:f} {:f} </lower_left>'\
                      '\n'.format(box[0][0], box[0][1])
        mesh_lines += '  <upper_right>{:f} {:f}</upper_right>'\
                      ''.format(box[1][0], box[1][1])
    else:
        mesh_lines += '  <dimension>{:d} {:d} {:d}</dimension>'\
                      '\n'.format(nmesh,nmesh,nmesh)
        mesh_lines += '  <lower_left>{:f} {:f} {:f}</lower_left>'\
                      '\n'.format(box[0][0], box[0][1], box[0][2])
        mesh_lines += '  <upper_right>{:f} {:f} {:f}</upper_right>'\
                      ''.format(box[1][0], box[1][1], box[1][2])
    replace_line('./tallies.xml', old, mesh_lines)

    return mesh, is2D

def get_adjoint(sp_file, mesh_dim):
    # Copy statepoint file if desired later
    shutil.copy(sp_file, sp_file[:-3] + '.adjoint.h5')
    FM = h5py.File(sp_file,'r') # Read the state point file
    # Normalize the tally result by the number of realizations
    FM_results = FM['tallies/tally 1/results'][()][:,:,0]/FM['tallies/tally 1/n_realizations'][()]
    FM_results_sd = FM['tallies/tally 1/results'][()][:,:,1]/FM['tallies/tally 1/n_realizations'][()]
    # Reshape the results obtain the fission matrix, 
    # mesh_dim is the number of bins in the mesh
    FM_results = FM_results.reshape([mesh_dim, mesh_dim]) 
    FM_results_sd = FM_results_sd.reshape([mesh_dim, mesh_dim]) 
    # get the number of neutrons produced in each mesh cell
    source = np.sum(FM_results,axis=0)
    source_sd = np.sum(FM_results_sd**2, axis=0)**(1/2)
    # The elements of the matrix need to be normalized by number of neutrons 
    # produced in each mesh cell
    FissionMatrix = copy.deepcopy(FM_results)
    print('source', source)
    for i in range(len(source)):
        if source[i] == 0:
            continue
        else:
            FissionMatrix[i,:] /= source[i]
    # Using the power method, obtain the dominant eigenvector of the fission 
    # matrix, this is the adjoint source
    adjoint = np.ones(mesh_dim)/mesh_dim # initialize the adjoint source
    diff = []
    n_iters = 300
    for i in range(n_iters):
        x1 = FissionMatrix@adjoint
        x1 = x1/np.sqrt(np.sum(x1**2))
        diff.append(np.sum((x1-adjoint)**2))
    adjoint = x1

    source_dict = {'source': source, 'source_sd': source_sd}
    
    with open('adjoint.pkl', 'wb') as f:
        pickle.dump(adjoint, f, pickle.HIGHEST_PROTOCOL)
    with open('source.pkl', 'wb') as f:
        pickle.dump(source_dict, f, pickle.HIGHEST_PROTOCOL)

    return adjoint

def make_sensitivity_tally(variable, nuclide, adjoint_filter, reaction=None, 
                           e_grid=None):
    sensitivity = openmc.Sensitivity() 
    sensitivity.variable = variable
    sensitivity.nuclide = nuclide
    if variable=='cross_section':
        if reaction is None:
            raise ValueError('Must specify reaction type for cross section\
                              sensitivity')
        sensitivity.reaction = reaction
    if e_grid is not None:
        sensitivity.energy = e_grid
    sens_tally = openmc.SensitivityTally() 
    sens_tally.filters = [adjoint_filter]
    sens_tally.sensitivity = sensitivity
    return sens_tally


def initialize_sensitivity_run(adjoint, mesh, nuclides, reactions, e_grid,
                               nbatches, nparticles):
    # Set number of particles and batches for sensitivity run
    settings = openmc.Settings.from_xml('settings.xml')
    settings.batches = nbatches
    settings.particles = nparticles
    settings.inactive = 20

    #particle tracking for debugging
    # settings.verbosity=10
    # batch = [73]
    # generations = [i for i in range(10)]
    # particles = [i for i in range(100)]
    # particle_list = []
    # for b in batch:
    #     for g in generations:
    #         for p in particles:
    #             particle_list.append(b)
    #             particle_list.append(g)
    #             particle_list.append(p)
    # settings.track = particle_list

    settings.export_to_xml()


    # use same mesh that was used in calculating adjoint
    # Create an importance-weighted filter, this takes in "adjoint" which is
    # the adjoint source calculated in the previous step, and the 
    # importance_mesh
    adjoint_filter = openmc.ImportanceFilter(adjoint, mesh) 
    sens_list = []

    # Load key nuclides if necessary

    if nuclides=='KEY':
        nucs = get_nucs_from_mat_xml('materials.xml')
        nuclides = []
        for key_nuc in KEY_NUCLIDES:
            if key_nuc in nucs:
                nuclides.append(key_nuc)

    for nuclide in nuclides:
        # Multipole
        sens_tally = make_sensitivity_tally('multipole', nuclide, adjoint_filter)
        sens_list.append(sens_tally)

        # Curvefit
        sens_tally = make_sensitivity_tally('curve_fit', nuclide, adjoint_filter)
        sens_list.append(sens_tally)

        # All reactions
        for rxn in reactions:
            sens_tally = make_sensitivity_tally('cross_section', nuclide, 
                    adjoint_filter, 
                    reaction=rxn,
                    e_grid=e_grid)
            sens_list.append(sens_tally)

    # Creates the sensitivity tally. This tally takes in the importance filter
    # and the sensitivity object.
    tallies = openmc.Tallies(sens_list)
    tallies.export_to_xml()

def run_sensitivity_sequence(xml_dir, nmesh, adj_nbatches, adj_nparticles,
                             clutch_nbatches, clutch_nparticles,
                             sens_nuclides, sens_rxns, sens_e_grid,
                             sens_run_dir=None):
    """
    xml_dir: path to directory with clean xmls
    nmesh: number of segments in each dimension for adjoint mesh (n x n x n)
    nbatchs: number of batches to run
    nparticles: number of particles per batch
    sens_nuclides: one of "KEY" or list of nuclides
    """

    sp_file = 'statepoint.{:d}.h5'.format(adj_nbatches)
    basedir = os.getcwd()
    e_grid = np.geomspace(1e-5, 2e7, 1000)
    clean_dir = basedir + '/' + xml_dir
    print('Clean xml directory: {:s}'.format(clean_dir))
    if sens_run_dir is None:
        sens_dir = 'sensitivity_run/' 
    mkdir(sens_dir)
    print('Sensitivity run directory : {:s}'.format(sens_dir))
    os.chdir(sens_dir)

    print('Initializing adjoint run')
    mesh, is2D = initialize_adjoint_run(clean_dir, nmesh, adj_nbatches, 
                                        adj_nparticles)
    print('Executing adjoint run')
    openmc.run()

    if is2D:
        mesh_dim = nmesh**2
    else:
        mesh_dim = nmesh**3 

    print('Loading adjoint from {:s}'.format(sp_file))
    adjoint = get_adjoint(sp_file, mesh_dim)

    print('Initializing sensitivity run')
    initialize_sensitivity_run(adjoint, mesh, sens_nuclides, sens_rxns, 
                               sens_e_grid, clutch_nbatches, clutch_nparticles)
    print('Executing sensitivity run')
    # Determine if we are on the cluster
    if re.search("cluster", socket.gethostname()):
        openmc.run(threads=16)
    else:
        openmc.run(threads=4)

    os.chdir(basedir)
