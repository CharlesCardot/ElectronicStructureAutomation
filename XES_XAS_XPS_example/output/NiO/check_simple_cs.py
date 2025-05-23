import numpy as np
import os
import sys
sys.path.insert(1,"/home/ccardot3/Python_Code/CharlesFunctions/")
import CharlesFunctions as CF
import json
from pathlib import Path

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.core import Structure

from collections import Counter

import warnings
warnings.filterwarnings("ignore")

'''
The goal of this script is to check to see that there is only one wyckoff positon for the
transition metal, and that it is surrounded by only a single ligand type.

This automation pipeline is not set up to properly handle multiple ligand types or more than
one wyckoff position, and if these two conditions are not satisfied the entire pipeline
should abort, because human intervention is required to appropriately address this system
'''

path_to_file = Path(os.path.realpath(__file__))
automation_index = next((i for i, p in enumerate(path_to_file.parts) if p == 'automation'), None)
parent_dir = path_to_file.parents[len(path_to_file.parts) - automation_index - 3]
# parent_dir should be one directory deeper than the "automation" directory, such as the automation_scratch directory
NAME = os.getcwd().split("/")[-1]

cif = parent_dir / str("materials/%s/%s.cif" % (NAME,NAME))
structure = Structure.from_file(cif)

dft_dir = str(parent_dir) + str("/output/%s/DFT/" % (NAME))
with open(parent_dir / str("materials/%s/%s.inp" % (NAME,NAME)),"r") as f:
    material_input = json.load(f)

with open(Path(dft_dir) / "out.scflow", "r") as f:
    scf_lines = f.readlines()

# Extract unique atom sites (Wyckoff positions)
elements_index = scf_lines.index(" El. sort        nr      rmin     default_rmax   factor     rmax \n")
end_elements_index = elements_index
for key,line in enumerate(scf_lines[elements_index:]):
    if line == "\n":
        end_elements_index = elements_index + key
        break
elements = [x.strip().split()[0:2] for x in scf_lines[elements_index+1:end_elements_index]]

################
simple_cs = True
################

def process_wyckoff(structure):
    """
    Finds every unique element site (Wyckoff site)
    in the structure and returns some useful data about
    the wyckoff postions. This uses the 'equivalent_atoms'
    part of the symmetry dataset from pymatgen's
    SpacegroupAnalyzer.

    Returns
    - unique_elem_counts: A dictionary which maps elements to the
    total number of wyckoff positions in the structure for that element.
    Ex: {'Cr': 1, 'P': 1, 'O': 3'}
    - wyckoff_site_indices: A dictionary which maps a element-number
    pair which denotes a wyckoff site (ex: "S2", Sulfur 2) to an
    index in the structure.
    Ex: {'Cr1': 0, 'P1': 5, 'O1': 8, 'O2': 16, 'O3': 20}
    - wyckoff_weights: A dictionary which maps element-number pairs
    which denote a wyckoff site to their symmetry multiplicity.
    Ex: {'Cr1: 1, 'P2': 1, 'O1': 2, 'O2': 1, 'O3': 1}
    """
    wyckoff_site_indices = {}
    unique_elem_counts = {}

    analyzer = SpacegroupAnalyzer(structure)
    symmetry_dataset = analyzer.get_symmetry_dataset()
    symmetry_equiv_atoms = symmetry_dataset['equivalent_atoms']

    for idx in set(symmetry_equiv_atoms):
        element = structure[idx].specie.symbol
        if element in unique_elem_counts:
            unique_elem_counts[element] += 1
            wyckoff_site_indices[element].append(idx)
        else:
            unique_elem_counts.update({element: 1})
            wyckoff_site_indices.update({element: [idx]})

    wyckoff_weights = {}
    elements_counter = {}
    divisor = min([val for idx, val in Counter(symmetry_equiv_atoms).items()])
    for idx, val in Counter(symmetry_equiv_atoms).items():
        element = structure[idx].specie.symbol
        if element in elements_counter:
            elements_counter[element] += 1
        else:
            elements_counter.update({element: 1})
        element = element + str(elements_counter[element])
        wyckoff_weights.update({element: round(val/divisor, 6)})

    return unique_elem_counts, wyckoff_site_indices, wyckoff_weights

unique_elem_counts, wyckoff_site_indices, wyckoff_weights = process_wyckoff(structure)

# Check for only a single TM wyckoff
if unique_elem_counts[material_input['TM']] > 1:
    simple_cs = False

# Check for multiple ligand types
center_atom_index = wyckoff_site_indices[material_input['TM']][0]
coords = structure.cart_coords
voronoi = VoronoiNN(tol=0.5)
voronoi_neighbors = voronoi.get_nn_info(structure, center_atom_index)
ligands = set([n['site'].species for n in voronoi_neighbors])
if len(ligands) > 1:
    simple_cs = False

if simple_cs:
    print("continue")
else:
    print("abort")
