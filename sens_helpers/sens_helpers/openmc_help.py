import openmc

def get_nucs_from_mat_xml(mat_xml):
    materials = openmc.Materials.from_xml(mat_xml)
    nucs = []
    for mat in materials:
        for nuc in mat.nuclides:
            nucs.append(nuc.name)
    nucs = list(set(nucs))
    return nucs
