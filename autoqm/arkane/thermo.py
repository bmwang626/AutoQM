#!/usr/bin/env python3
# encoding: utf-8

import os
import shutil
import subprocess

from .input_template.input_template import ArkaneThermo
from .species import make_arkane_species_input_file


def run_arkane_thermo(
    mol_id,
    smi,
    row_index,
    df,
    model_chemistry,
    subinputs_dir,
    suboutputs_dir,
    scratch_dir,
    RMG_path,
):
    # current dir
    current_dir = os.getcwd()

    # make scratch dir
    scratch_dir = os.path.join(scratch_dir, mol_id)
    os.makedirs(scratch_dir)

    # move to scratch dir
    os.chdir(scratch_dir)

    # make arkane species and thermo input files
    xyz = df.loc[row_index, "dft_xyz"]
    freq_path = df.loc[row_index, "freq_path"]
    sp_path = df.loc[row_index, "sp_path"]
    make_arkane_species_thermo_input_files(
        mol_id, smi, xyz, freq_path, sp_path, model_chemistry, use_bond_corrections=True
    )

    # run arkane
    subprocess.run(f"python {RMG_path}/Arkane.py input.py", shell=True)

    # move thermo file to suboutputs dir
    thermo_file = os.path.join("RMG_libraries", "thermo.py")
    shutil.copyfile(thermo_file, os.path.join(suboutputs_dir, f"{mol_id}_thermo.py"))

    # remove dummy input file
    tmp_input_file = os.path.join(subinputs_dir, f"{mol_id}.tmp")
    os.remove(tmp_input_file)

    # move back to current dir
    os.chdir(current_dir)

    # remove scratch dir
    shutil.rmtree(scratch_dir)


def make_arkane_species_thermo_input_files(
    mol_id, smi, xyz, freq_path, sp_path, model_chemistry, use_bond_corrections=True
):

    arkane_species_input_path = make_arkane_species_input_file(
        smi,
        xyz,
        freq_path,
        sp_path,
        model_chemistry,
        use_bond_corrections=use_bond_corrections,
    )

    # make arkane thermo input file
    arkane_thermo_input_file = "input.py"
    arkane_thermo_settings = {
        "model_chemistry": model_chemistry,
        "use_bond_corrections": use_bond_corrections,
        "species_label": mol_id,
        "species_file": arkane_species_input_path,
        "species_smiles": smi,
        "save_path": arkane_thermo_input_file,
    }

    ArkaneThermo(arkane_thermo_settings).save()
