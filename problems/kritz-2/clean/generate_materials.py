import openmc

mats = openmc.Materials()

mat = openmc.Material(1)
mat.name = "fuel"
mat.set_density('g/cm3', 9.58)
mat.add_nuclide('U235',     1.38925E-03, 'wo')
mat.add_nuclide('U238',     8.66894E-01, 'wo')
mat.add_nuclide('Pu239',    1.20886E-02, 'wo')
mat.add_nuclide('Pu240',    1.03982E-03, 'wo')
mat.add_nuclide('Pu241',    5.90909E-05, 'wo')
mat.add_nuclide('Pu242',    4.01728E-06, 'wo')
mat.add_element('O',        1.18487E-01, 'wo')
mat.add_nuclide('Am241',   3.82582E-05, 'wo')
mats.append(mat)

mat = openmc.Material(2)
mat.name = "zircaloy"
mat.set_density('g/cm3', 6.506)
mat.add_element('Zr', 0.98225, 'wo')
mat.add_element('Sn', 0.015, 'wo')
mat.add_element('Fe', 0.00125, 'wo')
mat.add_element('Cr', 0.001, 'wo')
mat.add_element('N', 0.0005, 'wo')
mats.append(mat)

mat = openmc.Material(3)
mat.name = "borated water"
mat.set_density('g/cm3', 0.99793)
mat.add_element('H'   , 0.11190 , 'wo')
mat.add_element('O',    0.88810 , 'wo')
mat.add_nuclide('B10',  8.847e-7, 'wo')
mat.add_nuclide('B11',  3.915e-6, 'wo')
mat.add_s_alpha_beta('c_H_in_H2O')
mats.append(mat)

mat = openmc.Material(4)
mat.name = "SS"
mat.set_density('g/cm3', 7.9)
mat.add_element('Cr'   , 0.20 , 'wo')
mat.add_element('Fe'   , 0.70 , 'wo')
mat.add_element('Ni'   , 0.10 , 'wo')
mats.append(mat)

mat = openmc.Material(5)
mat.name = "vapor"
mat.set_density('g/cm3', 1.8461e-5)
mat.add_element('O', 0.6667, 'wo')
mat.add_element('H', 0.3333, 'wo')
mats.append(mat)

mats.export_to_xml()
