#!/usr/bin/env python3
# encoding: utf-8

import numpy as np

from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions, GetStereoisomerCount

from rmgpy.qm.qmdata import QMData
from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.element import get_element
from rmgpy.qm.symmetry import PointGroupCalculator
from rmgpy.molecule.resonance import generate_kekule_structure

def enumerate_bonds(rmg_mol):
    '''
    Modified from ARC
    '''
    mol_list = generate_kekule_structure(rmg_mol)
    if mol_list:
        return mol_list[0].enumerate_bonds()
    else:
        return rmg_mol.enumerate_bonds()

def get_symbols_and_coords(xyz):
    lines = xyz.splitlines()
    symbols = list()
    coords = list()
    for line in lines:
        if line.strip():
            symbol, x, y, z = line.split()
            symbols.append(symbol)
            coords.append([float(x), float(y), float(z)])
    return symbols, np.array(coords)

def determine_symmetry(symbols, coords):
    """
    Modified from ARC
    """
    atom_numbers = list()
    for symbol in symbols:
        atom_numbers.append(get_element(symbol).number)
    # Coords is an N x 3 numpy.ndarray of atomic coordinates in the same order as `atom_numbers`.
    unique_id = '0'  # Just some name that the SYMMETRY code gives to one of its jobs.
    scr_dir = os.path.join('/tmp', 'symmetry_scratch')  # Scratch directory that the SYMMETRY code writes its files in.
    if not os.path.exists(scr_dir):
        os.makedirs(scr_dir)
    symmetry = optical_isomers = 1
    qmdata = QMData(
        groundStateDegeneracy=1,  # Only needed to check if valid QMData.
        numberOfAtoms=len(atom_numbers),
        atomicNumbers=atom_numbers,
        atomCoords=(coords, 'angstrom'),
        energy=(0.0, 'kcal/mol')  # Dummy
    )
    symmetry_settings = type('', (), dict(symmetryPath='symmetry', scratchDirectory=scr_dir))()
    pgc = PointGroupCalculator(symmetry_settings, unique_id, qmdata)
    pg = pgc.calculate()
    if pg is not None:
        symmetry = pg.symmetry_number
        optical_isomers = 2 if pg.chiral else optical_isomers
    return symmetry, optical_isomers

def get_optical_isomers(rdkit_mol, smi):
    """
    Modified from Xiaorui's easy_rmg_model
    """
    num_chiral = None
    chiral = Chem.FindMolChiralCenters(rdkit_mol, force=True, includeUnassigned=True)
    if chiral:
        try:
            opts = StereoEnumerationOptions(tryEmbedding=True, unique=True)
            isomers = EnumerateStereoisomers(rdkit_mol, opts)
        except:
            ## Todo: embed may not work for TSs
            print(f'Warning: {smi} needs manual check')
        else:
            num_chiral = len([isomer for isomer in isomers])
            if num_chiral < 1:
                print(f'Warning: {smi} needs manual check')
    else:
        num_chiral = 1
    return num_chiral