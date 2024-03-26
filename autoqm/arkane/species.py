#!/usr/bin/env python3
# encoding: utf-8

from rdkit import Chem

from rdmc.mol import RDKitMol
from rdmc.external.rmg import from_rdkit_mol

from arkane.statmech import is_linear

from .utils import enumerate_bonds, get_symbols_and_coords, determine_symmetry, get_optical_isomers
from .input_template.input_template import ArkaneSpecies

def make_arkane_species_input_file(smi, xyz, freq_path, sp_path, model_chemistry, use_bond_corrections=True, arkane_species_input_path='species.py'):
    # make RDKit molecule
    rdmc_mol = RDKitMol.FromSmiles(smi)
    rdkit_mol = rdmc_mol._mol

    # make RMG molecule
    rmg_mol = from_rdkit_mol(rdkit_mol)

    # get atom dict
    atom_dict = rmg_mol.get_element_count()

    # get bond dict
    if use_bond_corrections:
        bond_dict = enumerate_bonds(rmg_mol)
    else:
        bond_dict = {}

    # get multiplicity
    multiplicity = rmg_mol.multiplicity

    # get charge
    charge = rmg_mol.get_net_charge()

    # separate xyz into symbols and coords
    symbols, coords = get_symbols_and_coords(xyz)

    # get linear
    linear = is_linear(coordinates=coords)

    # get external symmetry and optical isomers
    external_symmetry, optical_isomers = determine_symmetry(symbols, coords)
    new_optical_isomers = get_optical_isomers(rdkit_mol, smi)
    if new_optical_isomers is not None:
        if new_optical_isomers >= 1:
            optical_isomers = new_optical_isomers
    
    # make arkane species input file
    arkane_species_settings = {
        'model_chemistry': model_chemistry,
        'use_bond_corrections': use_bond_corrections,
        'atom_dict': atom_dict,
        'bond_dict': bond_dict,
        'multiplicity': multiplicity,
        'charge': charge,
        'linear': linear,
        'external_symmetry': external_symmetry,
        'optical_isomers': optical_isomers,
        'freq': freq_path,
        'sp': sp_path,
        'save_path': arkane_species_input_path,
    }
    
    ArkaneSpecies(arkane_species_settings).save()

    return arkane_species_input_path