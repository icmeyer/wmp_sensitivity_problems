import numpy as np
import matplotlib.pyplot as plt

import openmc

import shutil
import os
import glob

XSLIB = openmc.data.DataLibrary.from_xml(os.environ.get('OPENMC_CROSS_SECTIONS'))

def clean_up(files):
    for filename in files:
        if os.path.exists(filename):
            os.remove(filename)

def mkdir(directory):
    if not os.path.exists(directory):
        os.mkdir(directory)

def scrub_material(mat):
    # Remove isotopes from material not present in xs library
    nucs_removed = []
    # Something super weird here: for some reason iterating using 
    # `for nuc in mat.nuclides` will skip some nuclides! Fixed
    # by converting mat.nuclides to list maybe worth reporting as 
    # an OpenMC bug
    nuclides_list = list(mat.nuclides)
    for nuc in nuclides_list:
        if XSLIB.get_by_material(nuc.name) is None:
            mat.remove_nuclide(nuc.name)
            nucs_removed.append(nuc.name)
    print('Removing nuclides not in xs library: ', nucs_removed)
    return mat

def make_inputs_fuel_T(T):
    burn_up_xml = 'materials_4.000000_MWd.xml'
    clean_up(['settings.xml', 'geometry.xml', 'materials.xml'])
    ### MATERIALS ###
    materials_used = openmc.Materials()

    # Get material compositions from converted CASMO file
    casmo_materials = openmc.Materials.from_xml(burn_up_xml)
    casmo_mat_dict = {}
    for casmo_mat in casmo_materials:
        casmo_mat = scrub_material(casmo_mat)
        casmo_mat_dict[casmo_mat.name] = casmo_mat
        materials_used.append(casmo_mat)
    # Rings are are numbered from inner to outer 
    # averaged fuel: fuel1, fuel2, ...
    # gd1: gd1_1, gd1_2, ...
    # gd2: gd2_1, gd2_2, ...

    # Gd Pincell
    gd_weight_ratio = 0.03
    gd = openmc.Material(name='gd')
    gd.add_element('Gd', 1.0)

    uo2 = openmc.Material(name='uo2')
    uo2.add_element('U', 1.0, enrichment=3.0, enrichment_type='wo')
    uo2.add_element('O', 2.0)

    gd_fuel = openmc.Material.mix_materials([uo2, gd], 
                                            [1-gd_weight_ratio, gd_weight_ratio],
                                            'wo',
                                            name='gd_fuel')
    gd_fuel.set_density('g/cc', 10.0)
    gd_fuel.temperature = T
    materials_used.append(gd_fuel)

    uo2_fuel = openmc.Material(name='uo2_fuel')
    uo2_fuel.add_element('U', 1.0, enrichment=3.0, enrichment_type='wo')
    uo2_fuel.add_element('O', 2.0)
    uo2_fuel.set_density('g/cc', 10.0)
    uo2_fuel.temperature = T
    materials_used.append(uo2_fuel)

    # Zirc2
    zirc = openmc.Material(name='zirc2')
    zirc.set_density('g/cc', 6.55)
    zirc.temperature = 300 + 273.15
    zirc.add_element('O', 0.006796, percent_type='ao')
    zirc.add_element('Cr', 0.001743, percent_type='ao')
    zirc.add_element('Fe', 0.001623, percent_type='ao')
    zirc.add_element('Ni', 0.000772, percent_type='ao')
    zirc.add_element('Zr', 0.978381, percent_type='ao')
    zirc.add_element('Sn', 0.010686, percent_type='ao')
    materials_used.append(zirc)

    # Moderator
    moderator = openmc.Material(name='moderator')
    moderator.set_density('g/cc', 0.74)
    moderator.temperature = 286 + 273.15
    moderator.add_element('H',  2.0)
    moderator.add_element('O',  1.0)
    materials_used.append(moderator)

    materials_used.export_to_xml()

    ### GEOMETRY ###
    nrings = 8
    rings = []
    fuel_radius = 0.5
    for i in range(nrings):
        rings.append(openmc.ZCylinder(r=fuel_radius*(i+1)/8))

    clad_outer_radius = openmc.ZCylinder(r=0.6)

    def make_cell(name, region, fill):
        return openmc.Cell(name=name, region=region, fill=fill)

    # Fill in fuel regions
    fuel_regions = []
    inner_region = -rings[0]
    fuel_regions.append(inner_region)
    for i in range(nrings-1):
        fuel_regions.append( (+rings[i] & -rings[i+1]) )
    gd_fuel_cells1 = []
    gd_fuel_cells2 = []
    uo2_fuel_cells = []
    for i in range(nrings):
        gd1_ring_mat = casmo_mat_dict['gd1_{:d}'.format(i + 1)]
        gd2_ring_mat = casmo_mat_dict['gd2_{:d}'.format(i + 1)]
        uo2_ring_mat = casmo_mat_dict['fuel{:d}'.format(i + 1)]
        gd_fuel_cells1.append(make_cell('gd_fuel1_{:d}'.format(i+1), 
                                        fuel_regions[i], gd1_ring_mat))
        gd_fuel_cells2.append(make_cell('gd_fuel2_{:d}'.format(i+1), 
                                        fuel_regions[i], gd2_ring_mat))
        uo2_fuel_cells.append(make_cell('uo2_fuel{:d}'.format(i+1), 
                                        fuel_regions[i], uo2_ring_mat))

    clad_region = +rings[nrings-1] & -clad_outer_radius

    clad_gd1 = openmc.Cell(name='gd_clad1', fill=zirc, region = clad_region)
    clad_gd2 = openmc.Cell(name='gd_clad2', fill=zirc, region = clad_region)
    clad_uo2 = openmc.Cell(name='uo2_clad', fill=zirc, region = clad_region)

    pitch = 1.6
    pin_cell_box = openmc.rectangular_prism(width=pitch, height=pitch)

    mod_region = pin_cell_box & +clad_outer_radius
    modcell_gd1 = openmc.Cell(name='gd_moderator', fill = moderator, region=mod_region)
    modcell_gd2 = openmc.Cell(name='gd_moderator', fill = moderator, region=mod_region)
    modcell_uo2 = openmc.Cell(name='uo2_moderator', fill = moderator, region=mod_region)

    gd_cells1 = gd_fuel_cells1 + [clad_gd1, modcell_gd1]
    gd_cells2 = gd_fuel_cells2 + [clad_gd2, modcell_gd2]
    uo2_cells = uo2_fuel_cells + [clad_uo2, modcell_uo2]
    g1 = openmc.Universe(name='gd_pin1', cells=gd_cells1)
    g2 = openmc.Universe(name='gd_pin2', cells=gd_cells2)
    u = openmc.Universe(name='uo2_pin', cells=uo2_cells)

    assembly = openmc.RectLattice(name='assembly')
    assembly.pitch = (pitch, pitch)
    assembly.lower_left = (-2*pitch, -2*pitch)
    assembly.universes = [[u, u, u, u],
                          [u, u, g1, u],
                          [u, g2, u, u],
                          [u, u, u, u]]

    assembly_region = openmc.rectangular_prism(width=4*pitch, height=4*pitch, 
                                              origin=[0,0], 
                                              boundary_type='reflective')
    assembly_cell = openmc.Cell(name='assembly_cell', 
                                fill=assembly,
                                region=assembly_region)

    root_universe = openmc.Universe(name='root universe', 
                                    cells=[assembly_cell])
    geometry = openmc.Geometry(root_universe)

    # Export to "geometry.xml"
    geometry.export_to_xml()

    ### SETTINGS ###
    batches = 1000
    inactive = 10
    particles = 10000

    # Instantiate a Settings object
    settings = openmc.Settings()
    settings.batches = batches
    settings.inactive = inactive
    settings.particles = particles
    settings.output = {'tallies': True}

    settings.temperature['method'] = 'nearest'
    settings.temperature['tolerance'] = 1000
    settings.temperature['multipole'] = True

    # Create an initial uniform spatial source distribution over fissionable zones
    bounds = [-0.63, -0.63, -0.63, 0.63, 0.63, 0.63]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    settings.source = openmc.Source(space=uniform_dist)

    # Export to "settings.xml"
    settings.export_to_xml()

if __name__=='__main__':
    T = 900
    make_inputs_fuel_T(T)
