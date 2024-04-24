"""
This module computes rate coefficient from .csv file containing energies and frequencies
"""

import ast
import logging
import os

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from arkane.kinetics import KineticsJob
from rmgpy import constants
from rmgpy.kinetics.tunneling import Eckart
from rmgpy.molecule import Molecule
from rmgpy.molecule.element import get_element
from rmgpy.reaction import Reaction
from rmgpy.species import Species, TransitionState
from utils import (
    get_lot_and_freq_scale,
    get_rmg_conformer,
    parse_command_line_arguments,
    xyz_str_to_coords,
)

# get logger
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

    rmg_conformers = []

    for spc_label in ["ts", "r1", "r2", "p1", "p2"]:

        if spc_label == "ts":
            multiplicity = row["multiplicity"]
        else:
            multiplicity = (
                Molecule().from_smiles(f"{row[f'{spc_label}smi']}").multiplicity
            )

        if spc_label == "ts":
            xyz = row["std_xyz_str"]
        else:
            xyz = row[f"{spc_label}_matched_std_xyz_str"]
        atomic_numbers, coords = xyz_str_to_coords(xyz)
        mass = (
            sum(
                [
                    get_element(int(atomic_number)).mass
                    for atomic_number in atomic_numbers
                ]
            )
            / constants.Na,
            "kg",
        )

        frequencies = ast.literal_eval(row[f"{spc_label}_dft_frequencies"])

        if energy_level == "qgdlpnoccsd(t)f12d/ccpvtzf12":
            e_electronic = row[f"{spc_label}_dlpno_sp_hartree"]
        elif energy_level == "qgwb97xd/def2svp":
            e_electronic = row[f"{spc_label}_dft_hartreefock_energy_hartree"]
        else:
            raise ValueError(f"Energy level {energy_level} not recognized")

        if isinstance(e_electronic, str):
            e_electronic = float(e_electronic)

        e_electronic = e_electronic * 2625500  # hartree to J/mol

        try:
            rmg_conformer = get_rmg_conformer(
                label=spc_label,
                level_of_theory=level_of_theory,
                e_electronic=e_electronic,
                frequencies=frequencies,
                coords=coords,
                numbers=atomic_numbers,
                mass=mass,
                multiplicity=multiplicity,
                freq_scale=freq_scale,
                molecule=None,
                use_atom_corrections=use_atom_corrections,
                use_bond_corrections=use_bond_corrections,
                scr_dir=scr_dir,
            )
            rmg_conformers.append(rmg_conformer)
        except Exception as e:
            logger.error(f"Error in getting rmg conformer for {row}: {e}")
            return None

    return rmg_conformers


def calc_rate_coefficient(
    row,
    freq_scale,
    energy_level,
    freq_level,
    energy_software,
    freq_software,
    Tmin=None,
    Tmax=None,
    Tcount=0,
    Tlist=None,
    sensitivity_conditions=None,
    three_params=True,
    scr_dir=None,
):

    output = get_rmg_conformer_from_df(
        row,
        freq_scale,
        energy_level,
        freq_level,
        energy_software=energy_software,
        freq_software=freq_software,
        use_atom_corrections=True,
        use_bond_corrections=False,
        scr_dir=scr_dir,
    )

    if output is None:
        return None

    ts, r1, r2, p1, p2 = output

    spc_r1 = Species().from_smiles(row["r1smi"])
    spc_r1.conformer = r1
    spc_r2 = Species().from_smiles(row["r2smi"])
    spc_r2.conformer = r2
    spc_p1 = Species().from_smiles(row["p1smi"])
    spc_p1.conformer = p1
    spc_p2 = Species().from_smiles(row["p2smi"])
    spc_p2.conformer = p2

    neg_frequency = row["neg_freq"]
    neg_frequency = (neg_frequency, "cm^-1")

    spc_ts = TransitionState(
        conformer=ts,
        frequency=neg_frequency,
        tunneling=Eckart(frequency=None, E0_reac=None, E0_TS=None, E0_prod=None),
    )

    rxn = Reaction(
        reactants=[spc_r1, spc_r2], products=[spc_p1, spc_p2], transition_state=spc_ts
    )
    kinetics_job = KineticsJob(
        reaction=rxn,
        Tmin=Tmin,
        Tmax=Tmax,
        Tcount=Tcount,
        Tlist=Tlist,
        sensitivity_conditions=sensitivity_conditions,
        three_params=three_params,
    )

    try:
        kinetics_job.generate_kinetics()
        return kinetics_job.reaction
    except Exception as e:
        logger.error(f"Error in generating kinetics for {row}: {e}")
        return None


def main():
    """
    The executable function
    """
    # 0. Parse input
    args = parse_command_line_arguments()

    df = pd.read_csv(args.csv_path)
    df = df.dropna()

    energy_level = args.energy_level
    freq_level = args.freq_level
    energy_software = args.energy_software
    freq_software = args.freq_software
    freq_scale = args.freq_scale

    reactions_list = Parallel(n_jobs=args.n_jobs, backend="multiprocessing")(
        delayed(calc_rate_coefficient)(
            row,
            freq_scale,
            energy_level,
            freq_level,
            energy_software,
            freq_software,
            scr_dir=f"./scratch_rxn_{idx}",
        )
        for idx, row in df.iterrows()
    )

    columns = ["A", "n", "Ea", "deltaHrxn", "deltaGrxn"]

    def get_values(rxn, key):
        if rxn is None:
            return None
        
        if key == "A":
            return rxn.kinetics.A.value_si
        elif key == "n":
            return rxn.kinetics.n.value_si
        elif key == "Ea":
            return rxn.kinetics.Ea.value_si
        elif key == "deltaHrxn":
            return rxn.get_enthalpy_of_reaction(298)
        elif key == "deltaGrxn":
            return rxn.get_free_energy_of_reaction(298)
        else:
            raise ValueError(f"Key {key} is not recognized")

    df["reaction"] = reactions_list

    for column in columns:
        df[column] = df["reaction"].apply(lambda x: get_values(x, column))

    df.to_csv(args.save_path, index=False)


if __name__ == "__main__":
    main()
