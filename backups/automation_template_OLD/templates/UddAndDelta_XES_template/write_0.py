import os
import sys
import json

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import Structure
from pathlib import Path

path_to_file = Path(os.path.realpath(__file__))
automation_index = next((i for i, p in enumerate(path_to_file.parts) if p == 'automation'), None)
parent_dir = path_to_file.parents[len(path_to_file.parts) - automation_index - 3]
# parent_dir should be one directory deeper than the "automation" directory, such as the automation_scratch directory

NAME = os.getcwd().split("/")[-1]
cif = parent_dir / str("materials/%s/%s.cif" % (NAME,NAME))
structure = Structure.from_file(cif)
sym_struct = SpacegroupAnalyzer(structure).get_symmetrized_structure()

# Changing number of k-points, dont need as many 
# for large unit cells
k_points = "@k@ 40 40 40"
if 100 <= sym_struct.lattice.volume < 200:
    k_points = "@k@ 30 30 30"
elif 200 <= sym_struct.lattice.volume < 300:
    k_points = "@k@ 20 20 20"
elif 300 <= sym_struct.lattice.volume:
    k_points = "@k@ 10 10 10"

def to_s(x):
    return f"{x:0.8f}"


filename = parent_dir / str("output/%s/0_Run%sFPLO.Quanty" % (NAME,NAME))
unique_sites = len(sym_struct.equivalent_sites)
with open(filename, "r") as f:
    lines = f.readlines()
for key, line in enumerate(lines):
    line = line.replace("NAME",NAME)
    line = line.replace("SPACEGROUP",str(sym_struct.spacegroup.int_number))
    line = line.replace("WYCHKOFF_NUMBER",str(unique_sites))
    
    lat_temp = tuple([to_s(x) for x in sym_struct.lattice.abc])
    line = line.replace("LATTICECONSTANTS",str(lat_temp).replace("(","").replace(")","").replace(",","").replace("'",""))
    angle_temp = tuple([to_s(x) for x in sym_struct.lattice.angles])
    line = line.replace("AXISANGLES",str(angle_temp).replace("(","").replace(")","").replace(",","").replace("'",""))
    
    line = line.replace("K-POINTS",k_points.replace("(","").replace(")","").replace(",","").replace("'",""))
    lines[key] = line


wyckoff_positions = []
for i, sites in enumerate(sym_struct.equivalent_sites):
    site = sites[0]
    row = ["@"+str(i+1)+"@", "".join([x for x in site.species_string if x.isalpha()])]
    row.append("@")
    row.extend([to_s(j) for j in site.frac_coords])
    wyckoff_positions.append(" ".join(row) + "\n")

pos = lines.index("WYCHKOFF_POSITIONS\n")
lines[pos:pos+1] = wyckoff_positions

with open(filename,"w") as f:
    for line in lines:
        f.write(line)





