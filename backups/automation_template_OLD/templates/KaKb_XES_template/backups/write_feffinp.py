import os
import sys
import json
import pymatgen as pm

import datetime

import numpy as np

from pathlib import Path

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import Structure
import pymatgen.io.feff
from pymatgen.io.feff.sets import FEFFDictSet
from pymatgen.io.feff.inputs import Header


path_to_file = Path(os.path.realpath(__file__))
automation_index = next((i for i, p in enumerate(path_to_file.parts) if p == 'automation'), None)
parent_dir = path_to_file.parents[len(path_to_file.parts) - automation_index - 3]
# parent_dir should be one directory deeper than the "automation" directory, such as the automation_scratch directory

NAME = os.getcwd().split("/")[-2]

cif = parent_dir / str("materials/%s/%s.cif" % (NAME,NAME))
structure = Structure.from_file(cif)

with open(parent_dir / str("materials/%s/%s.inp" % (NAME,NAME)),"r") as f:
    material_input = json.load(f)


feff_dict_set = FEFFDictSet(absorbing_atom = material_input["TM"], 
                            structure = structure, 
                            radius = 5,
                            config_dict = {})

filename = "feff.inp"
with open(filename,"r") as f:
    lines = f.readlines()


title_lines = str(Header(structure, comment = "Generated on " + str(datetime.datetime.now()) + " by Charles Cardot"))

# Fill in all the placeholders
for key, line in enumerate(lines):

    line = line.replace("TITLE", title_lines)

    line = line.replace("POTENTIALS", str(feff_dict_set.potential))
    line = line.replace("ATOMS", str(feff_dict_set.atoms))

    
    lines[key] = line


with open("feff.inp","w") as f:
    for line in lines:
        f.write(line)



