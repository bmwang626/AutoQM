import os
import shutil
import subprocess
import numpy as np

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions, GetStereoisomerCount

from rdmc.external.rmg import from_rdkit_mol

from arkane.statmech import is_linear

from rmgpy.qm.qmdata import QMData
from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.element import get_element
from rmgpy.qm.symmetry import PointGroupCalculator
from rmgpy.molecule.resonance import generate_kekule_structure

from radical_workflow.arkane.arkane_input import ArkaneSpecies, ArkaneThermo

def run_arkane_thermo(mol_id, smi, xyz, freq_path, sp_path, model_chemistry, subinputs_dir, suboutputs_dir, scratch_dir, RMG_path):
    # current dir
    current_dir = os.getcwd()

    # make scratch dir
    scratch_dir = os.path.join(scratch_dir, mol_id)
    os.makedirs(scratch_dir)

    # move to scratch dir
    os.chdir(scratch_dir)

    # make arkane input file
    make_arkane_input_file(mol_id, smi, xyz, freq_path, sp_path, model_chemistry, use_bond_corrections=True)

    # run arkane
    subprocess.run(f'python {RMG_path}/Arkane.py input.py', shell=True)

    # move thermo file to suboutputs dir
    thermo_file = os.path.join("RMG_libraries", "thermo.py")
    shutil.copyfile(thermo_file, os.path.join(suboutputs_dir, f'{mol_id}_thermo.py'))

    # remove dummy input file
    tmp_input_file = os.path.join(subinputs_dir, f'{mol_id}.tmp')
    os.remove(tmp_input_file)

    # move back to current dir
    os.chdir(current_dir)

    # remove scratch dir
    shutil.rmtree(scratch_dir)

def make_arkane_input_file(mol_id, smi, xyz, freq_path, sp_path, model_chemistry, use_bond_corrections=True):
    # make RDKit molecule
    rdkit_mol = Chem.MolFromSmiles(smi)

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
    
    # arkane species input path
    arkane_species_input_path = 'species.py'

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

    # arkane thermo input path
    arkane_thermo_input_file = 'input.py'
    arkane_thermo_settings = {
        'model_chemistry': model_chemistry,
        'use_bond_corrections': use_bond_corrections,
        'species_label': mol_id,
        'species_file': arkane_species_input_path,
        'species_smiles': smi,
        'save_path': arkane_thermo_input_file,
    }

    ArkaneThermo(arkane_thermo_settings).save()

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