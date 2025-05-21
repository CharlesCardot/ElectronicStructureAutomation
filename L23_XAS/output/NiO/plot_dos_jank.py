import numpy as np
import matplotlib.pyplot as plt
import os
import sys

from pathlib import Path

path_to_file = Path(os.path.realpath(__file__))
automation_index = next((i for i, p in enumerate(path_to_file.parts) if p == 'automation'), None)
parent_dir = path_to_file.parents[len(path_to_file.parts) - automation_index - 3]

utils_path = parent_dir.parents[0] / "utils"
sys.path.append(str(utils_path))
import CharlesFunctions as CF

NAME = "TbMn6Sn6"
dft_dir = "DFT"
material_input = {
    "TM": "Mn",
    "TM_VALENCE": "3d",
    "LIGAND": "Sn",
    "LIGAND_VALENCE": "4p"
}

print("NAME:", NAME)
print("dft_dir:", dft_dir)
print("material_input:\n", material_input)
print()

with open(Path(dft_dir) / "out.scf", "r") as f:
    scf_lines = f.readlines()

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
print()

TM_ldos = []
Ligand_ldos = []
for element in elements:
    for file in os.listdir(dft_dir):
        atom_num_pref = "00"
        if int(element[1]) > 9:
            atom_num_pref = "0"
        if int(element[1]) > 99:
            atom_num_pref = ""
        if file.startswith("+dos.sort" + atom_num_pref + element[1]):
            with open(Path(dft_dir) / file,"r") as f:
                lines = f.readlines()
            if material_input["TM_VALENCE"] in lines[0] and element[0] == material_input["TM"]:
                TM_ldos.append(np.loadtxt(Path(dft_dir) / file,skiprows=1))
            if material_input["LIGAND_VALENCE"] in lines[0] and element[0] == material_input["LIGAND"]:
                Ligand_ldos.append(np.loadtxt(Path(dft_dir) / file,skiprows=1))



TM_x = np.asarray(TM_ldos[0]).T[0]
TM_y = 0 
for entry in TM_ldos:
    TM_y += np.asarray(entry).T[1]
TM = np.asarray([TM_x, TM_y])
Ligand_x = np.asarray(Ligand_ldos[0]).T[0]
Ligand_y = 0 
for entry in Ligand_ldos:
    Ligand_y += np.asarray(entry).T[1]
Ligand_x = np.asarray(Ligand_ldos[0]).T[0]
Ligand = np.asarray([Ligand_x, Ligand_y])

TM = CF.Normalize(TM)
Ligand = CF.Normalize(Ligand)

plt.plot(TM[0], TM[1], label="TM")
plt.plot(Ligand[0], Ligand[1], label="Ligand")
plt.legend()
plt.show()


