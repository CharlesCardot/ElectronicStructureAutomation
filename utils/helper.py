"""
This file contains a lot of helper functions which
may be useful for high level setup up batch jobs or evaluating
which have completed.

It should NOT be used for any actual computation, plotting,
or complex file management. That functionality should be
contained within the template helper python files (ex: write_0.py).
"""

import json
import numpy as np
from pathlib import Path

from mp_api.client import MPRester
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import Composition
from pymatgen.core.periodic_table import Element


def get_mp_structures(mpids):
    """
    Materials Project and Pymatgen is trash, so this requires
    a dedicated conda env with the most updated versions of 
    pymatgen and mp-api to work.
    """


    MP_API_KEY = json.load(open("secrets.json","r"))["MP_API_KEY"]
    with MPRester(MP_API_KEY) as mpr:
        docs = mpr.materials.search(material_ids = mpids, fields = ['structure'])

    structures = [x.structure for x in docs]
    return structures



def generate_materials_inps(structures, output_path):
   
    for structure in structures:
        composition = Composition(str(structure.reduced_formula))
        tms = []
        for element in composition.elements:
            if element.is_transition_metal and element.Z < 31:
                tms.append(element.symbol)

        if len(tms) > 1:
            print("WARNING: multiple transition metals in structure, picking first one")
        tm = tms[0]

        for key, site in enumerate(structure.sites):
            if site.specie.symbol == tm:
                tm_index = key
                break

        nd = round(composition.oxi_state_guesses()[0][tm])

        coord_num = 0
        unique_ligands = set()
        r = 0.1
        ligand_count = 0
        while len(unique_ligands) < 2 and r < 10:
            neighbors = [x.specie.symbol for x in structure.get_neighbors(structure[tm_index], r)]
            for ligand in neighbors:
                if ligand in unique_ligands:
                    pass
                else:
                    unique_ligands.add(ligand)
                    if len(unique_ligands) == 1:
                        ligand_count = len(neighbors)

            r += 0.1
        
        unique_ligands.remove(tm)
        ligand = list(unique_ligands)[0]

        tm_valence = [x for x in Element(tm).electronic_structure.split(".") if "3d" in x][0][0:2]
        ligand_valence = Element(ligand).electronic_structure.split(".")[-1][0:2]
            
        input_file = {
             "NAME": str(structure.reduced_formula),
             "nd": nd,
             "coordination_num": ligand_count,
             "TM": tm,
             "TM_VALENCE": tm_valence,
             "LIGAND": ligand,
             "LIGAND_VALENCE": ligand_valence
        }
        print(input_file)

def setup_materials(structure, output_path):
    """
    Builds the materials folder in the main wrkdir of a automation
    job. This has the format of
    materials
        |-TiO2
            |- TiO2.cif
            |- TiO2.inp
        |-NiO
            |- NiO.cif
            |- NiO.inp
        |- ..
    
    """
    pass


generate_materials_inps(get_mp_structures(['mp-19006', 'mp-19009', 'mp-10695']), None)


