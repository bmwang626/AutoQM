import logging
import sys
sys.path.insert(0, "/home/gridsan/hwpang/Software/RMG-Py/")
import shutil
import os
import yaml
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from arkane.input import process_model_chemistry
from arkane.common import symbol_by_number, get_principal_moments_of_inertia
from arkane.encorr.corr import assign_frequency_scale_factor, get_atom_correction, get_bac
from arkane.ess import ess_factory
from arkane.modelchem import LOT, LevelOfTheory, CompositeLevelOfTheory, model_chem_to_lot
from arkane.thermo import ThermoJob

from rmgpy import constants
from rmgpy.kinetics.tunneling import Eckart
from rmgpy.molecule.element import get_element
from rmgpy.molecule.molecule import Molecule
from rmgpy.qm.qmdata import QMData
from rmgpy.qm.symmetry import PointGroupCalculator
from rmgpy.species import Species, TransitionState
from rmgpy.statmech import (IdealGasTranslation,
                            NonlinearRotor,
                            LinearRotor,
                            HarmonicOscillator,
                            Conformer)

def get_symmetry(coords, atom_numbers, scr_dir=None):

    scr_dir = scr_dir or os.path.join('.', 'scratch')
    os.makedirs(scr_dir, exist_ok=True)
    
    symmetry = optical_isomers = 1
    try:
        qmdata = QMData(
            groundStateDegeneracy=1,  # Only needed to check if valid QMData
            numberOfAtoms=len(atom_numbers),
            atomicNumbers=atom_numbers,
            atomCoords=(coords, 'angstrom'),
            energy=(0.0, 'kcal/mol')  # Only needed to avoid error
        )
        settings = type('', (),
                        dict(symmetryPath='symmetry',
                             scratchDirectory=scr_dir))()
        pgc = PointGroupCalculator(settings, '0', qmdata)  # '0' is an unique id used for calculator
        pg = pgc.calculate()
        if pg is not None:
            symmetry = pg.symmetry_number
            optical_isomers = 2 if pg.chiral else 1
            logging.debug(f"Symmetry algorithm found {optical_isomers} optical isomers "
                          f"and a symmetry number of {symmetry}")
        else:
            logging.warning("Symmetry algorithm errored when computing point group. "
                            "Using symmetry number=1 and optical isomers = 1 for "
                            "further calculations, which may not be true.")
        return symmetry, optical_isomers
    finally:
        shutil.rmtree(scr_dir)

def get_lot_and_freq_scale(energy_level: str,
                           freq_level: str,
                           energy_software: str = 'gaussian',
                           freq_software: str = 'gaussian',
                           freq_scale: Optional[float] = None):

    # Assign level of theory and frequency scale factor
    energy_method, energy_basis = energy_level.split('/')
    energy_lot = LevelOfTheory(method=energy_method, basis=energy_basis, software=energy_software)

    if energy_level == freq_level:
        level_of_theory = energy_lot
    else:
        freq_method, freq_basis = freq_level.split('/')
        freq_lot = LevelOfTheory(method=freq_method, basis=freq_basis, software=freq_software)
        level_of_theory = CompositeLevelOfTheory(freq=freq_lot, energy=energy_lot)
    return level_of_theory

def get_rotational_mode(coords, number, external_symmetry=None):
    if not external_symmetry:
        external_symmetry, _ = get_symmetry(coords, number)
    
    # Rotational
    moments_of_inertia = get_principal_moments_of_inertia(coords=coords,
                                                          numbers=number,)[0]
    if any([moment_of_inertia == 0.0 for moment_of_inertia in moments_of_inertia]):
        # this is a linear rotor
        moments_of_inertia = [moment_of_inertia for moment_of_inertia in moments_of_inertia
                              if moment_of_inertia != 0.0]
        if abs(moments_of_inertia[0] - moments_of_inertia[1]) > 0.01:
            raise Exceptions(f'Expected two identical moments of inertia for a linear rigis rotor, '
                             f'but got {moments_of_inertia}')
        return LinearRotor(inertia=(moments_of_inertia[0], "amu*angstrom^2"),
                           symmetry=external_symmetry)
    else:
        # this is a non-linear rotor
        return NonlinearRotor(inertia=(moments_of_inertia, "amu*angstrom^2"),
                              symmetry=external_symmetry)
    
def get_element_counts(number):
    # Get atoms count
    atoms = {}
    for atom_num in number:
        try:
            symbol = symbol_by_number[atom_num]
        except KeyError:
            raise ElementError('Could not recognize element number {0}.'.format(atom_num))
        atoms[symbol] = atoms.get(symbol, 0) + 1
    return atoms

def get_rmg_conformer_from_logs(label,
                                 energy_log_path,
                                 freq_log_path,
                                 level_of_theory,
                                 freq_scale=1,
                                 multiplicity=1,
                                 molecule=None,
                                 use_atom_corrections=True,
                                 use_bond_corrections=False,):
    
    energy_log = ess_factory(energy_log_path)
    freq_log = ess_factory(freq_log_path)
    
    # Get the coords, atom_numbers, and mass
    coords, number, mass = freq_log.load_geometry()
    # Get the symmetry info of the molecule
    external_symmetry, optical_isomers = get_symmetry(coords, number)
    
    # This script will automatically assign modes but with wrong E0 and unscaled freqs
    conformer, unscaled_frequencies = freq_log.load_conformer(symmetry=external_symmetry,
                                                              spin_multiplicity=multiplicity,
                                                              optical_isomers=optical_isomers,
                                                              label=label)
    
    conformer.coordinates = (coords, "angstroms")
    conformer.number = number
    conformer.mass = (mass, "amu")
    
    zpe_scale_factor = freq_scale / 1.014
    e_electronic = energy_log.load_energy(zpe_scale_factor)

    # Atom energy corrections
    if use_atom_corrections:
        # Get atoms count
        atoms = get_element_counts(number)
        atom_corrections = get_atom_correction(level_of_theory, atoms)
    else:
        atom_corrections = 0

    # Bond energy corrections
    if use_bond_corrections:
        # Get bonds count
        try:
            bonds = molecule.enumerate_bonds()
            bond_corrections = get_bac(level_of_theory, bonds, coords, number,
                                       bac_type='p', multiplicity=multiplicity)
        except AttributeError:
            raise ValueError('Cannot get BAC, since argument ``molecule`` is not provided.')
    else:
        bond_corrections = 0

    e_electronic_with_corrections = e_electronic + atom_corrections + bond_corrections
    zpe = freq_log.load_zero_point_energy() * zpe_scale_factor if len(number) > 1 else 0
    conformer.E0 = ((e_electronic_with_corrections + zpe) * 0.001, 'kJ/mol')

    # Correct the frequencies
    for mode in conformer.modes:
        if isinstance(mode, HarmonicOscillator):
            mode.frequencies = (np.array(unscaled_frequencies) * freq_scale, "cm^-1")
    
    return conformer, zpe

def get_rmg_conformer(label,
                      level_of_theory,
                      e_electronic,
                      frequencies,
                      coords,
                      numbers,
                      mass,
                      multiplicity=1,
                      freq_scale=1,
                      molecule=None,
                      use_atom_corrections=True,
                      use_bond_corrections=False,):
    
    external_symmetry, optical_isomers = get_symmetry(coords, numbers,)
    
    modes = []
    # Translational
    translation = IdealGasTranslation(mass=mass)
    modes.append(translation)
    
    # Rotational
    rotation = get_rotational_mode(coords, numbers,
                                   external_symmetry=external_symmetry)
    modes.append(rotation)
    
    # Vibrational
    frequencies = np.array(frequencies)
    frequencies = frequencies[frequencies >= 0]
    vibration = HarmonicOscillator(frequencies=(frequencies * freq_scale, "cm^-1"))
    modes.append(vibration)
    
    # Atom energy corrections
    if use_atom_corrections:
        atoms = get_element_counts(numbers)
        atom_corrections = get_atom_correction(level_of_theory, atoms)
    else:
        atom_corrections = 0

    # Bond energy corrections
    if use_bond_corrections:
        # Get bonds count
        try:
            bonds = molecule.enumerate_bonds()
            bond_corrections = get_bac(level_of_theory, bonds, coords, numbers,
                                       bac_type='p', multiplicity=multiplicity)
        except AttribureError:
            raise ValueError('Cannot get BAC, since argument ``molecule`` is not provided.')
    else:
        bond_corrections = 0
        
    e_electronic_with_corrections = e_electronic + atom_corrections + bond_corrections
    
    if len(numbers) > 1:
        zpe_scale_factor = freq_scale / 1.014
        scaled_zpe = 0.5 * constants.h * constants.c * constants.Na \
                     * np.sum(frequencies) * 100 * zpe_scale_factor
    else:
        scaled_zpe = 0

    e0 = ((e_electronic_with_corrections + scaled_zpe) * 0.001, 'kJ/mol')
    
    return Conformer(E0=e0, modes=modes,
                     spin_multiplicity=multiplicity,
                     optical_isomers=optical_isomers), scaled_zpe

# read input info

energy_log_path = '/home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/reactants_products_calculations/calculations/aug11b/output/DFT_opt_freq/outputs/outputs_0/id0.log'
freq_log_path = '/home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/reactants_products_calculations/calculations/aug11b/output/DFT_opt_freq/outputs/outputs_0/id0.log'
energy_level = 'QG-wb97xd/def2svp'
# energy_level = 'b3lyp/631g(d,p)'
energy_software = 'gaussian'
freq_level = 'QG-wb97xd/def2svp'
# freq_level = 'b3lyp/631g(d,p)'
freq_software = 'gaussian'
freq_scale = 0.986

# from input smiles csv
multiplicity = 1
smi = "[O:1]([O:2][H:4])[H:3]"
molecule = Molecule().from_smiles(smi)

level_of_theory = get_lot_and_freq_scale(energy_level=energy_level,
                                            freq_level=freq_level,
                                            energy_software=energy_software,
                                            freq_software=freq_software,
                                            freq_scale=freq_scale)

conf, zpe = get_rmg_conformer_from_logs('',
                            energy_log_path=energy_log_path,
                            freq_log_path=freq_log_path,
                            level_of_theory=level_of_theory,
                            freq_scale=freq_scale,
                            multiplicity=multiplicity,
                            molecule=molecule,
                            use_atom_corrections=True,
                            use_bond_corrections=True,)

spc = Species(molecule=[molecule])
spc.conformer = conf
thermo_job = ThermoJob(species=spc, thermo_class="NASA")
thermo_job.generate_thermo()
print(spc.thermo)

from rmgpy.kinetics.arrhenius import Arrhenius
from scipy.optimize import curve_fit

plt.figure()
Ts = np.arange(10, 800, 10)
Qs = np.array([spc.conformer.get_partition_function(T) for T in Ts])
arr = Arrhenius().fit_to_data(Ts, Qs, 's^-1')
print(arr)
plt.plot(Ts, Qs, label='Q')
plt.plot(Ts, [arr.get_rate_coefficient(T) for T in Ts], "--", label='arr')
plt.legend()
plt.yscale('log')
plt.savefig("Q.pdf")
# print(spc.thermo.to_thermo_data())
# rmg_conformer = get_rmg_conformer(label=hash_id,
#                                 level_of_theory=level_of_theory,
#                                 e_electronic=e_electronic,
#                                 frequencies=frequencies,
#                                 coords=coords,
#                                 numbers=atom_numbers,
#                                 mass=mass,
#                                 multiplicity=multiplicity,
#                                 freq_scale=freq_scale,
#                                 molecule=None,
#                                 use_atom_corrections=True,
#                                 use_bond_corrections=True,)


results = {
    "e0_zpe": conf.E0.value_si,
    "zpe": zpe,
    "h298": conf.get_enthalpy(298.15),
    "s298": conf.get_entropy(298.15),
    "g298": conf.get_free_energy(298.15),
    "q298": conf.get_partition_function(298.15),
    "cp298": conf.get_heat_capacity(298.15),
}

print(conf.E0.value_si)
print(conf.get_enthalpy(298.15))
print(conf.get_entropy(298.15))
print(conf.get_partition_function(298.15))
print(conf.get_heat_capacity(298.15))
