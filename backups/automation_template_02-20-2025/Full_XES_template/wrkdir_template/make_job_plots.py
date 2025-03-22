import sys
sys.path.insert(1,"/home/ccardot3/Python_Code/CharlesFunctions/")
import CharlesFunctions as CF

import numpy as np
import matplotlib.pyplot as plt
import scipy
import datetime
import re
import os

import plot_utils

from pathlib import Path
from matplotlib.ticker import ScalarFormatter

from pymatgen.vis.structure_vtk import StructureVis
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from PIL import Image

plot_params_dict = {
    "iso": {"color": "black"},
    "z_lin": {"color": "tab:orange"},
    "y_lin": {"color": "tab:green"},
    "x_lin": {"color": "tab:purple"},
    "z_rig": {"color": "blue"},
    "z_lef": {"color": "red"},
    "y_rig": {"color": "blue"},
    "y_lef": {"color": "red"},
    "x_rig": {"color": "blue"},
    "x_lef": {"color": "red"},
}

Ka_pp_dict = {
    "left": -30,
    "right": 5,
    "output_filename": "Ka_composite.png"
}

Kb_pp_dict = {
    "left": -30,
    "right": 5,
    "output_filename": "Kb_composite.png"
}

LDOS_pp_dict = {
    "sDOS": {"color": "tab:blue", "linewidth": 1.5},
    "pDOS": {"color": "tab:orange", "linewidth": 1.5},
    "dDOS": {"color": "tab:green", "linewidth": 1.5},
    "fDOS": {"color": "tab:red", "linewidth": 1.5},
    "totalDOS": {"color": "black", "linewidth": 1.5, "alpha": 0.6},
    "left": -100,
    "right": 100,
    "output_filename": "LDOS_composite.png"
}

VtC_pp_dict = {
    "z_dir_jd": {"color": "tab:red", "alpha": 0.5},
    "z_dir_jq": {"color": "tab:blue", "alpha": 0.5},
    "z_dir_dq": {"color": "black", "alpha": 0.5},
    "y_dir_jd": {"color": "tab:red", "alpha": 0.5},
    "y_dir_jq": {"color": "tab:blue", "alpha": 0.5},
    "y_dir_dq": {"color": "black", "alpha": 0.5},
    "x_dir_jd": {"color": "tab:red", "alpha": 0.5},
    "x_dir_jq": {"color": "tab:blue", "alpha": 0.5},
    "x_dir_dq": {"color": "black", "alpha": 0.5},
    "left": -50,
    "right": 50,
    "output_filename": "VtC_composite.png"
}

crystal_pp_dict = {
    "xdir_filename": "crystal_xdir.png",
    "ydir_filename": "crystal_ydir.png",
    "zdir_filename": "crystal_zdir.png",
}

BSandLDOS_pp_dict = {
    "emin": -10,
    "emax": 10,
    "output_filename": "bs_and_ldos.png"
}

## Material and File Management
# parent_dir should be one directory deeper than the "automation" directory, 
# such as the automation_scratch directory
path_to_file = Path(os.path.realpath(__file__))
automation_index = next((i for i, p in enumerate(path_to_file.parts) if p == 'automation'), None)
parent_dir = path_to_file.parents[len(path_to_file.parts) - automation_index - 3]

NAME = os.getcwd().split("/")[-1]
cif = parent_dir / str("materials/%s/%s.cif" % (NAME,NAME))
structure = Structure.from_file(cif)

import json
with open(parent_dir / str("materials/%s/%s.inp" % (NAME,NAME)),"r") as f:
    material_input = json.load(f)

# Add name to default plot_params
plot_params_dict.update({"NAME": NAME})

## DFT Due Diligence Plots
# SCF Convergence
plot_utils.make_SCFConvPlots("DFT/out.scflow", {"output_filename": "SCFconv_scflow.png"})
plot_utils.make_SCFConvPlots("DFT/out.scf", {"output_filename": "SCFconv_scf.png"})
image_paths = ["SCFconv_scflow.png", "SCFconv_scf.png"]
output_path = "stitched_SCFconv.png"
plot_utils.stitch_images_horizontal(image_paths, output_path)

# Band Structure, Wannier Bands, and LDOS
plot_utils.make_BandStructAndLdosPlots(Path("DFT"), material_input, BSandLDOS_pp_dict)

# Core to Core Spectra Plots
plot_utils.make_QuantyPlots("LFMcalc_Ka", {**Ka_pp_dict, **plot_params_dict})
plot_utils.make_QuantyPlots("LFMcalc_Kb", {**Kb_pp_dict, **plot_params_dict})

# FEFF Plots
with open(Path("VtC_FEFF/wyckoff_weights.json"), 'r') as f:
    wyckoff_weights = json.load(f)
VtC_pp_dict.update({"wyckoff_weights": wyckoff_weights})
LDOS_pp_dict.update({"wyckoff_weights": wyckoff_weights})

VtC_folders = [d for d in os.listdir(Path("VtC_FEFF")) if os.path.isdir(Path("VtC_FEFF") / d)]
VtC_folders.remove("VtC_template")
print("Found VtC Folders:", VtC_folders)

VtC_folders_elements = [re.sub(r'[0-9]', '', folder.split('_')[-1]) for folder in VtC_folders]
print(VtC_folders_elements)
VtC_folders_dict = {value: [index for index, element in enumerate(VtC_folders_elements) if element == value] for value in set(VtC_folders_elements)}
print(VtC_folders_dict)

# LDOS Plots
LDOS_avg_filenames = {}
for folder in VtC_folders:
    LDOS_pp_dict["output_filename"] = "LDOS" + folder.replace("VtC", "") + ".png"
    plot_utils.make_LDOSPlots([Path("VtC_FEFF") / folder], {**LDOS_pp_dict, **plot_params_dict})
for key, val in VtC_folders_dict.items():
    LDOS_pp_dict["output_filename"] = "LDOS_" + key + "_avg.png"
    LDOS_avg_filenames.update({key: LDOS_pp_dict["output_filename"]})
    plot_utils.make_LDOSPlots([Path("VtC_FEFF") / VtC_folders[i] for i in val], {**LDOS_pp_dict, **plot_params_dict})

# VtC Plots
VtC_avg_filenames = {}
for folder in VtC_folders:
    VtC_pp_dict["output_filename"] = folder + ".png"
    plot_utils.make_FEFFPlots([Path("VtC_FEFF") / folder], {**VtC_pp_dict, **plot_params_dict})
for key, val in VtC_folders_dict.items():
    VtC_pp_dict["output_filename"] = "VtC_" + key + "_avg.png"
    VtC_avg_filenames.update({key: VtC_pp_dict["output_filename"]})
    plot_utils.make_FEFFPlots([Path("VtC_FEFF") / VtC_folders[i] for i in val], {**VtC_pp_dict, **plot_params_dict})

# Crystal Visualization
plot_utils.make_crystalvisplots(structure, {"direction": "x", "output_filename": "crystal_xdir.png"})
plot_utils.make_crystalvisplots(structure, {"direction": "y", "output_filename": "crystal_ydir.png"})
plot_utils.make_crystalvisplots(structure, {"direction": "z", "output_filename": "crystal_zdir.png"})
image_paths = ["crystal_zdir.png", "crystal_ydir.png", "crystal_xdir.png"]
output_path = "stitched_crystal_image.png"
plot_utils.stitch_images_horizontal(image_paths, output_path)

# Meta Data
text = f"""
METADATA

Plots Generated by: Charles Cardot
Date: {str(datetime.datetime.now())[0:-5]}
System: {NAME}

NOTATION
(axes of propegation, axes of polarization, helicity)
[axes of propegation, dipole, quadrupole]

ex: (z,0,+1) = photon propegating along z-axis with
+1 (right) circular polarization.

ex: [x,d,q] = photon propegating along the x-axis
which is an average over all y-z polarizations, and
includes contributions from both dipole and quadrupole
transitions.
"""
plot_utils.text_to_image(text = text, output_path = "header.png", font_size=60)

image_paths = ["header.png", "Ka_composite.png", "Kb_composite.png"]
for key, val in VtC_avg_filenames.items():
    image_paths.append(val)
    image_paths.append(LDOS_avg_filenames[key])
image_paths = image_paths + ["stitched_crystal_image.png", "stitched_SCFconv.png", "bs_and_ldos.png"]
output_pdf_path = f"{NAME}_summary.pdf"
plot_utils.convert_images_to_pdf(image_paths, output_pdf_path)
plot_utils.move_png_files(folder_name = f"{NAME}_plots")


