#!/usr/bin/env python3
# encoding: utf-8

import os
import shutil
import subprocess

from .input_template.input_template import Arkanekinetics
from .species import make_arkane_species_input_file

def run_arkane_kinetics(ts_id, rxn_smi, row_index, df, model_chemistry, subinputs_dir, suboutputs_dir, scratch_dir, RMG_path):
    # current dir
    current_dir = os.getcwd()

    # make scratch dir
    scratch_dir = os.path.join(scratch_dir, ts_id)
    os.makedirs(scratch_dir)

    # move to scratch dir
    os.chdir(scratch_dir)

    # make arkane input file for reactants and products
    for spc in ['r1', 'r2', 'p1', 'p2']:
        spc_smi = df.loc[row_index, f'{spc}_smi']
        spc_xyz = df.loc[row_index, f'{spc}_dft_xyz']
        spc_freq_path = df.loc[row_index, f'{spc}_freq_path']
        spc_sp_path = df.loc[row_index, f'{spc}_sp_path']
        make_arkane_species_input_file(spc_smi, spc_xyz, spc_freq_path, spc_sp_path, model_chemistry, use_bond_corrections=False, arkane_species_input_path=f'{spc}.py')
    ts_xyz = df.loc[row_index, 'ts_dft_xyz']
    ts_freq_path = df.loc[row_index, 'ts_freq_path']
    ts_sp_path = df.loc[row_index, 'ts_sp_path']
    make_arkane_reaction_kinetics_input_files(rxn_smi, ts_xyz, ts_freq_path, ts_sp_path, model_chemistry)

    # run arkane
    subprocess.run(f'python {RMG_path}/Arkane.py input.py', shell=True)

    # move kinetics file to suboutputs dir
    kinetics_file = os.path.join("RMG_libraries", "reactions.py")
    shutil.copyfile(kinetics_file, os.path.join(suboutputs_dir, f'{ts_id}.py'))

    # remove dummy input file
    tmp_input_file = os.path.join(subinputs_dir, f'{ts_id}.tmp')
    os.remove(tmp_input_file)

    # move back to current dir
    os.chdir(current_dir)

    # remove scratch dir
    shutil.rmtree(scratch_dir)

def make_arkane_reaction_kinetics_input_files(rxn_smi, ts_xyz, ts_freq_path, ts_sp_path, model_chemistry):

    # make arkane input file for TS
    r_complex_smi, p_complex_smi = rxn_smi.split('>>')
    ts_smi = r_complex_smi # use reactant complex smi as TS smi
    arkane_ts_input_path = make_arkane_species_input_file(ts_smi, ts_xyz, ts_freq_path, ts_sp_path, model_chemistry, use_bond_corrections=False, arkane_species_input_path='TS.py')

    reactants = {f'r{i}': f'r{i}.py' for i in range(1, 3)}
    products = {f'p{i}': f'p{i}.py' for i in range(1, 3)}
    arkane_kinetics_input_file = 'input.py'
    arkane_kinetics_settings = {
        'model_chemistry': model_chemistry,
        'use_bond_corrections': False,
        'reactants': reactants,
        'products': products,
        'TS': arkane_ts_input_path,
        'tunneling': 'Eckart',
        'save_path': arkane_kinetics_input_file
    }

    Arkanekinetics(arkane_kinetics_settings).save()
