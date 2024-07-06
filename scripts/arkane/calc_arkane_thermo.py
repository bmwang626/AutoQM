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
    multiplicity,
    xyz_str,
    frequencies,
    energy,
    freq_scale,
    energy_level,
    freq_level,
    energy_software,
    freq_software,
    use_atom_corrections,
    use_bond_corrections,
    bac_type,
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

    atomic_numbers, coords = xyz_str_to_coords(xyz_str)
    mass = (
        sum([get_element(int(atomic_number)).mass for atomic_number in atomic_numbers])
        / constants.Na,
        "kg",
    )
    frequencies = ast.literal_eval(frequencies)
    frequencies = np.array(frequencies)
    e_electronic = energy * 2625500  # hartree to J/mol

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
            bac_type=bac_type,
            scr_dir=scr_dir,
        )
    except Exception as e:
        logger.error(f"Error in getting rmg conformer for multiplicity {multiplicity}, xyz_str {xyz_str}, frequencies {frequencies}, energy {energy}, freq_scale {freq_scale}, energy_level {energy_level}, freq_level {freq_level}, energy_software {energy_software}, freq_software {freq_software}, use_atom_corrections {use_atom_corrections}, use_bond_corrections {use_bond_corrections}, molecule {molecule}, scr_dir {scr_dir}: {e}")
        return None

    return rmg_conformer


def calc_thermo(
    smi,
    multiplicity,
    xyz_str,
    frequencies,
    energy,
    freq_scale,
    energy_level,
    freq_level,
    energy_software,
    freq_software,
    use_bond_corrections,
    scr_dir=None,
):

    molecule = Molecule().from_smiles(smi)

    rmg_conformer = get_rmg_conformer_from_df(
        multiplicity=multiplicity,
        xyz_str=xyz_str,
        frequencies=frequencies,
        energy=energy,
        freq_scale=freq_scale,
        energy_level=energy_level,
        freq_level=freq_level,
        energy_software=energy_software,
        freq_software=freq_software,
        use_atom_corrections=True,
        use_bond_corrections=use_bond_corrections,
        molecule=molecule,
        scr_dir=scr_dir,
        bac_type="p",
    )

    if rmg_conformer is None:
        return None, None

    spc = Species(molecule=[molecule])
    spc.conformer = rmg_conformer
    thermo_job = ThermoJob(species=spc, thermo_class="wilhoit")
    thermo_job.generate_thermo()
    p_thermo = thermo_job.species.thermo

    molecule = Molecule().from_smiles(smi)

    rmg_conformer = get_rmg_conformer_from_df(
        multiplicity=multiplicity,
        xyz_str=xyz_str,
        frequencies=frequencies,
        energy=energy,
        freq_scale=freq_scale,
        energy_level=energy_level,
        freq_level=freq_level,
        energy_software=energy_software,
        freq_software=freq_software,
        use_atom_corrections=True,
        use_bond_corrections=use_bond_corrections,
        molecule=molecule,
        scr_dir=scr_dir,
        bac_type="m",
    )

    if rmg_conformer is None:
        return None, None

    spc = Species(molecule=[molecule])
    spc.conformer = rmg_conformer
    thermo_job = ThermoJob(species=spc, thermo_class="wilhoit")
    thermo_job.generate_thermo()
    m_thermo = thermo_job.species.thermo

    return p_thermo, m_thermo


def main():
    """
    The executable function
    """
    # 0. Parse input
    args = parse_command_line_arguments()

    setattr(args, "freq_column", "species_dft_frequencies")

    if args.energy_level == "qgwb97xd/def2svp":
        energy_column = "species_dft_hartreefock_energy_hartree"
    elif args.energy_level == "qgdlpnoccsd(t)f12d/ccpvtzf12":
        energy_column = "species_dlpno_sp_hartree"
    else:
        raise ValueError(f"Energy level {args.energy_level} is not supported")

    df = pd.read_csv(args.csv_path)
    energy_level = args.energy_level
    freq_level = args.freq_level
    energy_software = args.energy_software
    freq_software = args.freq_software
    freq_scale = args.freq_scale

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

    thermoss = Parallel(n_jobs=args.n_jobs)(
        delayed(calc_thermo)(
            row["asmi"],
            row["multiplicity"],
            row["std_xyz_str"],
            row["species_dft_frequencies"],
            row[energy_column],
            freq_scale,
            energy_level,
            freq_level,
            energy_software,
            freq_software,
            use_bond_corrections=not args.no_bac_for_thermo,
            scr_dir=args.scratch_dir / f"thermo_{idx}",
        )
        for idx, row in tqdm(df.iterrows())
    )

    def get_values(thermo, key):
        if thermo is None:
            return None

        if key == "H298":
            return thermo.get_enthalpy(298)
        elif key == "S298":
            return thermo.get_entropy(298)
        elif key == "Cp300":
            return thermo.get_heat_capacity(300)
        elif key == "Cp400":
            return thermo.get_heat_capacity(400)
        elif key == "Cp500":
            return thermo.get_heat_capacity(500)
        elif key == "Cp600":
            return thermo.get_heat_capacity(600)
        elif key == "Cp800":
            return thermo.get_heat_capacity(800)
        elif key == "Cp1000":
            return thermo.get_heat_capacity(1000)
        elif key == "Cp1500":
            return thermo.get_heat_capacity(1500)
        else:
            raise ValueError(f"Key {key} is not supported")


    df["p_thermo"] = [thermos[0] for thermos in thermoss]
    df["m_thermo"] = [thermos[1] for thermos in thermoss]

    for column in columns:
        df["p_" + column] = df["p_thermo"].apply(lambda x: get_values(x, column))
        df["m_" + column] = df["m_thermo"].apply(lambda x: get_values(x, column))

    df.to_csv(args.save_path, index=False)


if __name__ == "__main__":
    main()
