import os
import sys
sys.path.insert(1,"/home/ccardot3/Python_Code/CharlesFunctions/")
import CharlesFunctions as CF
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

path_to_file = Path(os.path.realpath(__file__))
automation_index = next((i for i, p in enumerate(path_to_file.parts) if p == 'automation'), None)
parent_dir = path_to_file.parents[len(path_to_file.parts) - automation_index - 3]
# parent_dir should be one directory deeper than the "automation" directory, such as the automation_scratch directory

NAME = os.getcwd().split("/")[-1]
dft_dir = str(parent_dir) + str("/output/%s/DFT/" % (NAME))
with open(parent_dir / str("materials/%s/%s.inp" % (NAME,NAME)),"r") as f:
    material_input = json.load(f)

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

#plt.plot(TM[0], TM[1], label="TM")
#plt.plot(Ligand[0], Ligand[1], label="Ligand")
#plt.show()

'''
Find appropriate emin and emax

This is going to be super jank but ... 

- Start by seperating the total dos out into chunks
- Use the integral under the curve to determine how 'important' a chunk is
- Apply some cutoff to the chunks (ex: must have at least 10% weight and be within -15 < e < 15)
- Use that to find the largest emax and smallest emin
'''

chunk_weight_cutoff = 0.05
edge_buffer = 0.5 # Added to emax, subtracted from emin, to give a buffer
W = 2 # width of gaussian supression
left_window_edge = -8
right_window_edge = 8
min_chunk_height = 0.001

# Transition Metal (TM) Chunking
TM_chunks = []
inchunk = False
chunk = [[],[]]
orbital_center = 0
for key,val in enumerate(TM[0]):
    orbital_center += TM[1][key]/np.sum(TM[1]) * val
    if TM[1][key] > min_chunk_height:
        chunk[0].append(val)
        chunk[1].append(TM[1][key])
        inchunk = True
    elif inchunk:
        TM_chunks.append(np.asarray(chunk))
        chunk = [[],[]]
        inchunk = False
    if inchunk and key == len(TM[0]) - 1:
        TM_chunks.append(np.asarray(chunk))

 
# Weights chunks that are too far away from Fermi energy lower
# This is probably overely complicated
from scipy import integrate
import copy

TM_chunks_temp = copy.deepcopy(TM_chunks)
for chunk_key,chunk in enumerate(TM_chunks_temp):
    for key,val in enumerate(chunk[0]):
        if val < left_window_edge:
            chunk[1][key] = chunk[1][key] * CF.Gaussian(val, N=1, W=W, center=orbital_center + left_window_edge)
        elif val > right_window_edge:
            chunk[1][key] = chunk[1][key] * CF.Gaussian(val, N=1, W=W, center=orbital_center + right_window_edge)
    TM_chunks_temp[chunk_key] = chunk 
TM_chunks = TM_chunks_temp

print("\nTransition Metal Chunks")
print("TM Orbital Center", orbital_center)
for key,chunk in enumerate(TM_chunks):
    print("Chunk:", key, "Left Edge:", chunk[0][0], "Right Edge:", chunk[0][-1], "Weight:", integrate.trapz(chunk[1],chunk[0]))
    label = "chunk " + str(key)
    #plt.plot(chunk[0],chunk[1], label=label)

TM_chunks = [chunk for chunk in TM_chunks if (integrate.trapz(chunk[1],chunk[0]) > chunk_weight_cutoff)]
TM_emin = np.round(np.min(np.concatenate([chunk[0] for chunk in TM_chunks])) - edge_buffer,2)
TM_emax = np.round(np.max(np.concatenate([chunk[0] for chunk in TM_chunks])) + edge_buffer,2)
print("TM_emin", TM_emin)
print("TM_emax", TM_emax)


# Ligand Chunking
Ligand_chunks = []
inchunk = False
chunk = [[],[]]
orbital_center = 0
for key,val in enumerate(Ligand[0]):
    orbital_center += Ligand[1][key]/np.sum(Ligand[1]) * val
    if Ligand[1][key] > min_chunk_height:
        chunk[0].append(val)
        chunk[1].append(Ligand[1][key])
        inchunk = True
    elif inchunk:
        Ligand_chunks.append(np.asarray(chunk))
        chunk = [[],[]]
        inchunk = False
    if inchunk and key == len(Ligand[0]) - 1:
        Ligand_chunks.append(np.asarray(chunk))


# Weights chunks that are too far away from Fermi energy lower
# This is probably overely complicated
from scipy import integrate
import copy
Ligand_chunks_temp = copy.deepcopy(Ligand_chunks)
for chunk_key,chunk in enumerate(Ligand_chunks_temp):
    for key,val in enumerate(chunk[0]):
        if val < left_window_edge:
            chunk[1][key] = chunk[1][key] * CF.Gaussian(val, N=1, W=W, center=orbital_center + left_window_edge)
        elif val > right_window_edge:
            chunk[1][key] = chunk[1][key] * CF.Gaussian(val, N=1, W=W, center=orbital_center + right_window_edge)
    Ligand_chunks_temp[chunk_key] = chunk 
Ligand_chunks = Ligand_chunks_temp

print("\nLigand Chunks")
print("Ligand Orbital Center", orbital_center)
for key,chunk in enumerate(Ligand_chunks):
    print("Chunk:", key, "Left Edge:", chunk[0][0], "Right Edge:", chunk[0][-1], "Weight:", integrate.trapz(chunk[1],chunk[0]))
    label = "chunk " + str(key)
    #plt.plot(chunk[0],chunk[1], label=label)

Ligand_chunks = [chunk for chunk in Ligand_chunks if (integrate.trapz(chunk[1],chunk[0]) > chunk_weight_cutoff)]
Ligand_emin = np.round(np.min(np.concatenate([chunk[0] for chunk in Ligand_chunks])) - edge_buffer,2)
Ligand_emax = np.round(np.max(np.concatenate([chunk[0] for chunk in Ligand_chunks])) + edge_buffer,2)
print("Ligand_emin", Ligand_emin)
print("Ligand_emax", Ligand_emax)

#plt.legend()
#plt.show()

emin = str(min(TM_emin, Ligand_emin))
emax = str(max(TM_emax, Ligand_emax))


# Find all the atom sites that FPLO identifies
elements_index = scf_lines.index("No.  Element WPS CPA-Block    X                      Y                      Z\n")
end_elements_index = elements_index
for key,line in enumerate(scf_lines[elements_index:]):
    if line == "\n":
        end_elements_index = elements_index + key
        break
elements = [x.strip().split()[0:2] for x in scf_lines[elements_index+1:end_elements_index]]
print("\nAtom Sites")
print(elements)

TM_start = str(np.min([int(x[0]) for x in elements if material_input["TM"] in x]))
TM_end = str(np.max([int(x[0]) for x in elements if material_input["TM"] in x]))
Ligand_start = str(np.min([int(x[0]) for x in elements if material_input["LIGAND"] in x]))
Ligand_end = str(np.max([int(x[0]) for x in elements if material_input["LIGAND"] in x]))

# Writing everything out to WF.Quanty file
filename = parent_dir / str("output/%s/1_Run%sFPLO_WF.Quanty" % (NAME,NAME))
with open(filename,"r") as f:
    lines = f.readlines()
for key, line in enumerate(lines):
    line = line.replace("NAME",NAME)
    line = line.replace("EMAX",emax)
    line = line.replace("EMIN",emin)
    line = line.replace("TM_START",TM_start)
    line = line.replace("TM_END",TM_end)
    line = line.replace("TM_ELEMENT",material_input["TM"])
    line = line.replace("TM_VALENCE",material_input["TM_VALENCE"])
    line = line.replace("LIGAND_START",Ligand_start)
    line = line.replace("LIGAND_END",Ligand_end)
    line = line.replace("LIGAND_ELEMENT",material_input["LIGAND"])
    line = line.replace("LIGAND_VALENCE",material_input["LIGAND_VALENCE"])
    lines[key] = line

with open(filename,"w") as f:
    for line in lines:
        f.write(line)













