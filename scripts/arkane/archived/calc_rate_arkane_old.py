"""
This module computes rate coefficient from .csv file containing energies and frequencies
"""

import argparse
import ast
import os

import numpy as np
import pandas as pd
from rmgpy import constants
from rmgpy.kinetics import Arrhenius
from rmgpy.kinetics.tunneling import Eckart
from rmgpy.molecule.element import get_element
from rmgpy.reaction import Reaction
from utils import (
    get_lot_and_freq_scale,
    get_rmg_conformer,
    parse_command_line_arguments,
)


def get_rmg_conformer_from_df(
    row,
    freq_scale,
    energy_level,
    freq_level,
    energy_software,
    freq_software,
    use_atom_corrections,
    use_bond_corrections,
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
            prefix = ""
        else:
            prefix = spc_label + "_"

        multiplicity = row[f"{prefix}multiplicity"]
        xyz = ast.literal_eval(row[f"{prefix}std_xyz"])
        xyz = np.array(xyz)
        coords = xyz[:, 3:]
        atom_numbers = xyz[:, 1]
        mass = (
            sum([get_element(int(atom_number)).mass for atom_number in atom_numbers])
            / constants.Na,
            "kg",
        )
        frequencies = ast.literal_eval(row[f"{prefix}frequencies"])
        frequencies = np.array(frequencies)
        e_electronic = row[f"{prefix}hf"] * 2625500  # hartree to J/mol

        rmg_conformer = get_rmg_conformer(
            label=spc_label,
            level_of_theory=level_of_theory,
            e_electronic=e_electronic,
            frequencies=frequencies,
            coords=coords,
            numbers=atom_numbers,
            mass=mass,
            multiplicity=multiplicity,
            freq_scale=freq_scale,
            molecule=None,
            use_atom_corrections=use_atom_corrections,
            use_bond_corrections=use_bond_corrections,
        )
        rmg_conformers.append(rmg_conformer)

    return rmg_conformers


def get_E0min_and_Q(Ts, conformers, option="ms"):
    if option not in ["min", "max", "ms"]:
        raise

    if option == "min":
        E0_min = min([conf.E0.value_si for conf in conformers])
        min_idx = np.argmin([conf.E0.value_si for conf in conformers])
        E0_ref = E0_min
        Q = np.zeros_like(Ts)
        for idx in range(Ts.shape[0]):
            T = Ts[idx]
            Q[idx] = np.array([conformers[min_idx].get_partition_function(T)])

    elif option == "max":
        E0_max = max([conf.E0.value_si for conf in conformers])
        max_idx = np.argmax([conf.E0.value_si for conf in conformers])
        E0_ref = E0_max
        Q = np.zeros_like(Ts)
        for idx in range(Ts.shape[0]):
            T = Ts[idx]
            Q[idx] = np.array([conformers[max_idx].get_partition_function(T)])

    elif option == "ms":
        E0_min = min([conf.E0.value_si for conf in conformers])
        E0_ref = E0_min
        Uis = np.array([conf.E0.value_si - E0_ref for conf in conformers])
        Q = np.zeros_like(Ts)
        for idx in range(Ts.shape[0]):
            T = Ts[idx]
            qis = np.array([conf.get_partition_function(T) for conf in conformers])
            Q[idx] = np.sum(qis * np.exp(-Uis / constants.R / T))

    return E0_ref, Q


def calc_rate_coefficient(
    row,
    Temps,
    freq_scale,
    energy_level,
    freq_level,
    energy_software,
    freq_software,
):

    ts, r1, r2, p1, p2 = get_rmg_conformer_from_df(
        row,
        freq_scale,
        energy_level,
        freq_level,
        energy_software=energy_software,
        freq_software=freq_software,
        use_atom_corrections=True,
        use_bond_corrections=False,
    )

    rxn = Reaction(reactants=[r1, r2], products=[p1, p2], transitionState=ts)
    job = KineticsJob(
        reaction=rxn,
        Tmin=Tmin,
        Tmax=Tmax,
        Tcount=Tcount,
        Tlist=Tlist,
        sensitivity_conditions=sensitivity_conditions,
        three_params=three_params,
    )
    # def calc_rate_coefficient(reactants, products):
    #     E0_reac = 0
    #     Q_reac = 1
    #     for reactant in reactants:
    #         E0_min, Q = get_E0min_and_Q(Temps, [reactant], option='min')
    #         E0_reac += E0_min
    #         Q_reac *= Q

    #     E0_prod = 0
    #     for product in products:
    #         E0_min, _ = get_E0min_and_Q(Temps, [product], option='min')
    #         E0_prod += E0_min

    #     E0_TS, Q_TS = get_E0min_and_Q(Temps, [ts], option='min')

    #     neg_frequency = (ts.frequency.value_si, 'cm^-1')

    #     eckart = Eckart(frequency = neg_frequency,
    #                     E0_reac=(E0_reac, 'J/mol'),
    #                     E0_TS=(E0_TS, 'J/mol'),
    #                     E0_prod=(E0_prod, 'J/mol'))

    #     dE0 = E0_TS - E0_reac

    #     ks_w_tunnel = np.zeros_like(Temps)
    #     ks_wo_tunnel = np.zeros_like(Temps)

    #     for idx in range(Temps.shape[0]):
    #         T = Temps[idx]
    #         ks_wo_tunnel[idx] = (constants.kB * T / constants.h * Q_TS[idx] / Q_reac[idx]) * np.exp(-dE0 / constants.R / T)
    #         ks_w_tunnel[idx] = ks_wo_tunnel[idx] * eckart.calculate_tunneling_factor(T)
    #     return (ks_w_tunnel, ks_wo_tunnel)

    # fwd_ks_w_tunnel, fwd_ks_wo_tunnel = calc_rate_coefficient([r1, r2], [p1, p2])
    # rev_ks_w_tunnel, rev_ks_wo_tunnel = calc_rate_coefficient([p1, p2], [r1, r2])

    return fwd_ks_w_tunnel, fwd_ks_wo_tunnel, rev_ks_w_tunnel, rev_ks_wo_tunnel


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

    save_path = os.path.join(os.path.dirname(args.csv_path), "rate_coefficients.csv")

    fwd_rxn_order = 2
    rev_rxn_order = 2
    k_units_dict = {1: "s^-1", 2: "m^3/(mol*s)", 3: "m^6/(mol^2*s)"}

    Temps = np.linspace(298, 3000, 200)

    df[
        [
            "fwd_w_tunnel_A",
            "fwd_w_tunnel_Ea",
            "fwd_wo_tunnel_A",
            "fwd_wo_tunnel_Ea",
            "rev_w_tunnel_A",
            "rev_w_tunnel_Ea",
            "rev_wo_tunnel_A",
            "rev_wo_tunnel_Ea",
        ]
    ] = None

    for idx, row in df.iterrows():

        fwd_ks_w_tunnel, fwd_ks_wo_tunnel, rev_ks_w_tunnel, rev_ks_wo_tunnel = (
            calc_rate_coefficient(
                row,
                Temps,
                freq_scale,
                energy_level,
                freq_level,
                energy_software,
                freq_software,
            )
        )

        fwd_ks_w_tunnel_fit = Arrhenius().fit_to_data(
            Temps,
            fwd_ks_w_tunnel,
            kunits=k_units_dict[fwd_rxn_order],
            three_params=False,
        )
        fwd_ks_wo_tunnel_fit = Arrhenius().fit_to_data(
            Temps,
            fwd_ks_wo_tunnel,
            kunits=k_units_dict[fwd_rxn_order],
            three_params=False,
        )
        rev_ks_w_tunnel_fit = Arrhenius().fit_to_data(
            Temps,
            rev_ks_w_tunnel,
            kunits=k_units_dict[rev_rxn_order],
            three_params=False,
        )
        rev_ks_wo_tunnel_fit = Arrhenius().fit_to_data(
            Temps,
            rev_ks_wo_tunnel,
            kunits=k_units_dict[rev_rxn_order],
            three_params=False,
        )

        fwd_ks_w_tunnel_A = fwd_ks_w_tunnel_fit.A.value_si
        fwd_ks_w_tunnel_Ea = fwd_ks_w_tunnel_fit.Ea.value_si
        fwd_ks_wo_tunnel_A = fwd_ks_wo_tunnel_fit.A.value_si
        fwd_ks_wo_tunnel_Ea = fwd_ks_wo_tunnel_fit.Ea.value_si
        rev_ks_w_tunnel_A = rev_ks_w_tunnel_fit.A.value_si
        rev_ks_w_tunnel_Ea = rev_ks_w_tunnel_fit.Ea.value_si
        rev_ks_wo_tunnel_A = rev_ks_wo_tunnel_fit.A.value_si
        rev_ks_wo_tunnel_Ea = rev_ks_wo_tunnel_fit.Ea.value_si

        df.loc[
            idx,
            [
                "fwd_w_tunnel_A",
                "fwd_w_tunnel_Ea",
                "fwd_wo_tunnel_A",
                "fwd_wo_tunnel_Ea",
                "rev_w_tunnel_A",
                "rev_w_tunnel_Ea",
                "rev_wo_tunnel_A",
                "rev_wo_tunnel_Ea",
            ],
        ] = [
            fwd_ks_w_tunnel_A,
            fwd_ks_w_tunnel_Ea,
            fwd_ks_wo_tunnel_A,
            fwd_ks_wo_tunnel_Ea,
            rev_ks_w_tunnel_A,
            rev_ks_w_tunnel_Ea,
            rev_ks_wo_tunnel_A,
            rev_ks_wo_tunnel_Ea,
        ]

    df.to_csv(save_path, index=False)


if __name__ == "__main__":
    main()
