import os
import re
import sys
import json
import shutil

import datetime

import numpy as np

from pathlib import Path
from collections import Counter

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import Structure
import pymatgen.io.feff
from pymatgen.io.feff.sets import FEFFDictSet
from pymatgen.io.feff.inputs import Header

"""
This script is designed to identify every wyckoff position in a crystal
and then create a new VtC-XES feff_template.inp for that atom at that wyckoff
position. It will also create and label VtC folders acording to the 
element-wyckoff number (ex: VtC_Fe1) for each wyckoff position.
"""

######
# parent_dir should be one directory deeper than the "automation" directory, 
# such as the automation_scratch directory
#####
path_to_file = Path(os.path.realpath(__file__))
automation_index = next((i for i, p in enumerate(path_to_file.parts) if p == 'automation'), None)
parent_dir = path_to_file.parents[len(path_to_file.parts) - automation_index - 3]

NAME = os.getcwd().split("/")[-2]
cif = parent_dir / str("materials/%s/%s.cif" % (NAME,NAME))
structure = Structure.from_file(cif)

# Do this just in case there is only 1 absorbing atom in the unit cell
structure = structure * np.asarray([2,2,2])

with open(parent_dir / str("materials/%s/%s.inp" % (NAME,NAME)),"r") as f:
    material_input = json.load(f)

sg_analyzer = SpacegroupAnalyzer(structure)
symmetry_operations = sg_analyzer.get_symmetry_operations()

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
print("Unique Elements:", unique_elem_counts)
print("Wyckoff Site Indices:", wyckoff_site_indices)
print("Wyckoff Weights:", wyckoff_weights)

with open("wyckoff_weights.json", "w") as f:
    json.dump(wyckoff_weights, f, indent=4)

# Have to do this for FEFFDictSet to work with element_index
for site in structure:
    site.species = site.specie.symbol

for element, unique_num in unique_elem_counts.items():
    for i in range(unique_num):
        element_index = int(wyckoff_site_indices[element][i])
        feff_dict_set = FEFFDictSet(absorbing_atom = element_index, 
                            structure = structure, 
                            radius = 7,
                            config_dict = {})

        filename = "feff_template.inp"
        with open(filename,"r") as f:
            lines = f.readlines()
  
        title_lines = str(Header(structure, comment = "Generated on " + str(datetime.datetime.now()) + " by Charles Cardot"))
        title_lines = "\n".join([line for key, line in enumerate(title_lines.split("\n")) if line[0] != "*" or key == 0])
        
        # Fill in all the placeholders
        for key, line in enumerate(lines):
        
            line = line.replace("TITLE", title_lines)
        
            line = line.replace("POTENTIALS", str(feff_dict_set.potential))
            line = line.replace("ATOMS", str(feff_dict_set.atoms))
        
            
            lines[key] = line
        
        # Create VtC folder for this element-number pair
        dest_folder = "VtC_" + element + str(i + 1)
        shutil.copytree("VtC_template", dest_folder)
        
        with open(Path(dest_folder) / "feff_template.inp","w") as f:
            for line in lines:
                f.write(line)

