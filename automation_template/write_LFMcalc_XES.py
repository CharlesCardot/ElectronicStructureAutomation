import os
import sys
import json
import pymatgen as pm

import datetime

import numpy as np
import pandas as pd

from pathlib import Path

path_to_file = Path(os.path.realpath(__file__))
automation_index = next((i for i, p in enumerate(path_to_file.parts) if p == 'automation'), None)
parent_dir = path_to_file.parents[len(path_to_file.parts) - automation_index - 3]
# parent_dir should be one directory deeper than the "automation" directory, such as the automation_scratch directory

NAME = os.getcwd().split("/")[-2]
DFT_PATH = "/".join(os.getcwd().split("/")[0:-1] + ["DFT"])

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import Structure

cif = parent_dir / str("materials/%s/%s.cif" % (NAME,NAME))
structure = Structure.from_file(cif)

with open(Path(DFT_PATH) / "out.scf", "r") as f:
    scf_lines = f.readlines()

with open(parent_dir / str("materials/%s/%s.inp" % (NAME,NAME)),"r") as f:
    material_input = json.load(f)

# Find the quanty file to edit
quanty_files = [filename for filename in os.listdir(os.getcwd()) if filename.endswith(".Quanty") or filename.endswith(".quanty")]
if len(quanty_files) > 1:
    print("Multiple quanty files detected, which one to edit?")
    print("Exiting ...")
    quit()
filename = quanty_files[0]

with open(filename,"r") as f:
    lines = f.readlines()

# Extract unique atom sites (Wyckoff positions)
elements_index = scf_lines.index(" El. sort        nr      rmin     default_rmax   factor     rmax \n")
end_elements_index = elements_index
for key,line in enumerate(scf_lines[elements_index:]):
    if line == "\n":
        end_elements_index = elements_index + key
        break
elements = [x.strip().split()[0:2] for x in scf_lines[elements_index+1:end_elements_index]]
print("Wyckoff positions")
print(elements)

# Find the FPLO number for the TM Wyckoff position
element_symbols = [x[0] for x in elements]
# Need to update workflow to include materials with 
# more than one unique Wyckoff position for the Transition Metal
if element_symbols.count(material_input["TM"]) > 1:
    print("More than one unique TM Wyckoff position, aborting calculation")
    quit()

# Find the Wyckoff number for one of the TM atoms in the calculation
for element in elements:
    atom_num_pref = "00"
    if int(element[1]) > 9:
        atom_num_pref = "0"
    if int(element[1]) > 99:
        atom_num_pref = ""
    if element[0] == material_input["TM"]:
        TM_NUMBER = atom_num_pref + element[1]
        break


# Extract general atom sites in the unit cell that FPLO identifies
elements_index = scf_lines.index("No.  Element WPS CPA-Block    X                      Y                      Z\n")
end_elements_index = elements_index
for key,line in enumerate(scf_lines[elements_index:]):
    if line == "\n":
        end_elements_index = elements_index + key
        break
elements = [x.strip().split()[0:2] + x.strip().split()[4:7] for x in scf_lines[elements_index+1:end_elements_index]]
print("\nAtom Sites")
print(elements)

TMs = [x for x in elements if material_input["TM"] in x]
TM_FPLO_XYZ = " , ".join(TMs[0][2:5])

print("\nTM_FPLO_XYZ")
print(TM_FPLO_XYZ)
coordination_num = material_input["coordination_num"]
for bond_dist in np.arange(0.1,10,0.1):
    # Convert to angstroms
    arr = structure.get_sites_in_sphere(pt=np.asarray([float(x) * 0.529 for x in TMs[0][2:5]]), r=bond_dist)
    if len(arr) == coordination_num + 1:
        # Conver back to Bohr
        cluster_radius = str(np.round(2.1 * bond_dist,2))
        break
print("\ncluster_radius") 
print(cluster_radius)


zeta_df = pd.read_csv(parent_dir.parents[0] / "utils/All3d_zeta_values_tape.csv")

if int(material_input["nd"]) < 1:
    configuration = "2p63p64s03d1" # Technically, no 3d spin-orbit coupling if no 3d electrons
else:
    configuration = "2p63p64s03d" + str(material_input["nd"])
filtered_row = zeta_df[(zeta_df['Element'] == material_input["TM"]) & (zeta_df['Configuration'] == configuration)]
TM_ZETA_3D = str(filtered_row["3d"].iloc[0])

configuration = "2p53p64s03d" + str(material_input["nd"])
filtered_row = zeta_df[(zeta_df['Element'] == material_input["TM"]) & (zeta_df['Configuration'] == configuration)]
TM_ZETA_2P = str(filtered_row["2p"].iloc[0])

configuration = "2p63p54s03d" + str(material_input["nd"])
filtered_row = zeta_df[(zeta_df['Element'] == material_input["TM"]) & (zeta_df['Configuration'] == configuration)]
TM_ZETA_3P = str(filtered_row["3p"].iloc[0])

with open(parent_dir.parents[0] / str("utils/TM_LIFETIMES.json"), "r") as f:
    lifetimes = json.load(f)

# Fill in all the placeholders
for key, line in enumerate(lines):

    line = line.replace("-- Material:", "-- Material: " + NAME)
    line = line.replace("-- Date:", "-- Date: " + str(datetime.datetime.now()))

    line = line.replace("DFT_PATH",DFT_PATH)
    line = line.replace("TM_ELEMENT",material_input["TM"])
    line = line.replace("TM_NUMBER",TM_NUMBER)
    line = line.replace("TM_FPLO_XYZ",TM_FPLO_XYZ)
    line = line.replace("NUMBER_OF_D",str(material_input["nd"]))
    line = line.replace("TM_CLUSTER_RADIUS",cluster_radius)

    line = line.replace("TM_ZETA_2P",TM_ZETA_2P)
    line = line.replace("TM_ZETA_3P",TM_ZETA_3P)
    line = line.replace("TM_ZETA_3D",TM_ZETA_3D)

    lines[key] = line


with open(filename,"w") as f:
    for line in lines:
        f.write(line)



