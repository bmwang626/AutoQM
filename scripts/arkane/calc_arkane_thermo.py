"""
This module computes rate coefficient from .csv file containing energies and frequencies
"""

import ast
import os
import logging
from tqdm import tqdm

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from arkane.thermo import ThermoJob
from rmgpy import constants
from rmgpy.molecule.element import get_element
from rmgpy.molecule.molecule import Molecule
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

from utils import (
    get_lot_and_freq_scale,
    get_rmg_conformer,
    parse_command_line_arguments,
    xyz_str_to_coords,
)

logger = logging.getLogger()


def get_rmg_conformer_from_df(
    row,
    freq_scale,
    energy_level,
    freq_level,
    energy_software,
    freq_software,
    use_atom_corrections,
    use_bond_corrections,
    molecule,
    scr_dir=None,
):
    # No bond correction at this moment

    level_of_theory, freq_scale = get_lot_and_freq_scale(
        energy_level=energy_level,
        freq_level=freq_level,
        energy_software=energy_software,
        freq_software=freq_software,
        freq_scale=freq_scale,
    )

    multiplicity = row[f"multiplicity"]
    atomic_numbers, coords = xyz_str_to_coords(row[f"std_xyz_str"])
    mass = (
        sum([get_element(int(atomic_number)).mass for atomic_number in atomic_numbers])
        / constants.Na,
        "kg",
    )
    frequencies = ast.literal_eval(row[f"species_dft_frequencies"])
    frequencies = np.array(frequencies)
    e_electronic = row[f"species_dlpno_sp_hartree"] * 2625500  # hartree to J/mol

    try:
        rmg_conformer = get_rmg_conformer(
            label="",
            level_of_theory=level_of_theory,
            e_electronic=e_electronic,
            frequencies=frequencies,
            coords=coords,
            numbers=atomic_numbers,
            mass=mass,
            multiplicity=multiplicity,
            freq_scale=freq_scale,
            molecule=molecule,
            use_atom_corrections=use_atom_corrections,
            use_bond_corrections=use_bond_corrections,
            scr_dir=scr_dir,
        )
    except Exception as e:
        logger.error(f"Error in getting rmg conformer for {row}: {e}")
        return None

    return rmg_conformer


def calc_thermo(
    row,
    freq_scale,
    energy_level,
    freq_level,
    energy_software,
    freq_software,
    scr_dir=None,
):
    smi = row["asmi"]

    molecule = Molecule().from_smiles(smi)

    rmg_conformer = get_rmg_conformer_from_df(
        row=row,
        freq_scale=freq_scale,
        energy_level=energy_level,
        freq_level=freq_level,
        energy_software=energy_software,
        freq_software=freq_software,
        use_atom_corrections=True,
        use_bond_corrections=True,
        molecule=molecule,
        scr_dir=scr_dir,
    )

    if rmg_conformer is None:
        return None

    spc = Species(molecule=[molecule])
    spc.conformer = rmg_conformer
    thermo_job = ThermoJob(species=spc, thermo_class="NASA")
    thermo_job.generate_thermo()

    return thermo_job.species.thermo


def main():
    """
    The executable function
    """
    # 0. Parse input
    args = parse_command_line_arguments()

    df = pd.read_csv(args.csv_path)
    energy_level = args.energy_level
    freq_level = args.freq_level
    energy_software = args.energy_software
    freq_software = args.freq_software
    freq_scale = args.freq_scale

    save_path = os.path.join(
        os.path.dirname(args.csv_path),
        os.path.basename(args.csv_path).split(".csv")[0] + "_thermos.csv",
    )

    columns = [
        "H298",
        "S298",
        "Cp300",
        "Cp400",
        "Cp500",
        "Cp600",
        "Cp800",
        "Cp1000",
        "Cp1500",
    ]
    df[columns] = None

    thermos = Parallel(n_jobs=args.n_jobs)(
        delayed(calc_thermo)(
            row,
            freq_scale,
            energy_level,
            freq_level,
            energy_software,
            freq_software,
            scr_dir=f"./scratch_thermo_{idx}",
        )
        for idx, row in tqdm(df.iterrows())
    )

    for idx, thermo in tqdm(enumerate(thermos)):
        if thermo is None:
            continue

        df.loc[idx, columns] = [
            thermo.get_enthalpy(298),
            thermo.get_entropy(298),
            thermo.get_heat_capacity(300),
            thermo.get_heat_capacity(400),
            thermo.get_heat_capacity(500),
            thermo.get_heat_capacity(600),
            thermo.get_heat_capacity(800),
            thermo.get_heat_capacity(1000),
            thermo.get_heat_capacity(1500),
        ]

    df.to_csv(save_path, index=False)


if __name__ == "__main__":
    main()
