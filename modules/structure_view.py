import numpy as np
import pandas as pd
import plotly.graph_objs as go
from dash_bio.utils import PdbParser, create_mol3d_style


def make_structure_info(structure_id, AF3_data):
    parser = PdbParser(f"https://raw.githubusercontent.com/slt666666/CS_public_pdb/refs/heads/main/pdb/{str(structure_id).upper()}.PDB")
    structure_data = parser.mol3d_data()
    structure_styles = create_mol3d_style(
        structure_data['atoms'], visualization_type='cartoon', color_element='residue'
    )
    
    tmp_AF3_data = AF3_data[AF3_data["Gene ID"] == str(structure_id).upper()]
    length = tmp_AF3_data["Amino acids"].values[0]
    color_per_residue = np.repeat("purple", length+1)
    CC = tmp_AF3_data.CC.values[0].split("-")
    color_per_residue[int(CC[0])-1:int(CC[1])-1] = "gray"
    NB = tmp_AF3_data.NB.values[0].split("-")
    color_per_residue[int(NB[0])-1:int(NB[1])-1] = "orange"
    LRR = tmp_AF3_data.LRR.values[0].split("-")
    color_per_residue[int(LRR[0])-1:int(LRR[1])-1] = "blue"
    color_per_residue = list(color_per_residue)

    for i, (style, atom) in enumerate(zip(structure_styles, structure_data['atoms'])):
        # Note: we iterate over all atoms but have one color per residue.
        # We therefore get the residue index for every atom.
        residue_idx = atom['residue_index']

        style['color'] = color_per_residue[residue_idx]
    
    return structure_styles, structure_data, tmp_AF3_data.Remarks.values[0]
