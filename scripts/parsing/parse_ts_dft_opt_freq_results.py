#!/usr/bin/env python
# coding: utf-8
import os
import sys
import pandas as pd
import pickle as pkl
from tqdm import tqdm
from joblib import Parallel, delayed
from argparse import ArgumentParser

from autoqm.parser.dft_opt_freq_parser import dft_opt_freq_parser

parser = ArgumentParser()
parser.add_argument(
    "--input_smiles_path",
    type=str,
    required=True,
    help="path to a .csv file containing input smiles and ids",
)
parser.add_argument(
    "--output_file_name",
    type=str,
    required=True,
    help="name of the output file",
)
parser.add_argument(
    "--n_jobs",
    type=int,
    required=True,
    help="number of jobs to run in parallel",
)
args = parser.parse_args()

input_smiles_path = args.input_smiles_path
output_file_name = args.output_file_name
n_jobs = args.n_jobs

# input_smiles_path = "reactants_products_wb97xd_and_xtb_opted_ts_combo_results_hashed_chart_aug11b.csv"
# n_jobs = 8

df = pd.read_csv(input_smiles_path)
rxn_ids = df["id"].tolist()

log_paths = []
for rxn_id in rxn_ids:
    ids = rxn_id // 1000
    log_path = os.path.join(
        "output",
        "DFT_opt_freq",
        "outputs",
        f"rxns_{ids}",
        f"rxn_{rxn_id}",
        f"{rxn_id}.log",
    )
    log_paths.append(log_path)

out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(
    delayed(dft_opt_freq_parser)(path, is_ts=True, check_connectivity=False,) for path in tqdm(log_paths) # not able to use check_connectivity=True for TS
)

failed_jobs = dict()
valid_jobs = dict()
for rxn_id, (failed_job, valid_job) in zip(rxn_ids, out):
    if failed_job:
        failed_jobs[rxn_id] = failed_job
    if valid_job:
        valid_jobs[rxn_id] = valid_job

with open(os.path.join(f"{output_file_name}.pkl"), "wb") as outfile:
    pkl.dump(valid_jobs, outfile, protocol=pkl.HIGHEST_PROTOCOL)

with open(os.path.join(f"{output_file_name}_failed.pkl"), "wb") as outfile:
    pkl.dump(failed_jobs, outfile, protocol=pkl.HIGHEST_PROTOCOL)

print(f"Total number of molecules: {len(rxn_ids)}")
print(f"Total number of failed molecules: {len(failed_jobs)}")
print(failed_jobs)

rxn_id_to_DFT_opted_xyz_std_ori = {}
for rxn_id, valid_job in valid_jobs.items():
    rxn_id_to_DFT_opted_xyz_std_ori[rxn_id] = valid_job["dft_xyz_std_ori"]

with open(os.path.join(f"{output_file_name}_xyz_std_ori.pkl"), "wb") as outfile:
    pkl.dump(rxn_id_to_DFT_opted_xyz_std_ori, outfile, protocol=pkl.HIGHEST_PROTOCOL)

rxn_id_to_DFT_opted_xyz_input_ori = {}
for rxn_id, valid_job in valid_jobs.items():
    rxn_id_to_DFT_opted_xyz_input_ori[rxn_id] = valid_job["dft_xyz_input_ori"]

with open(os.path.join(f"{output_file_name}_xyz_input_ori.pkl"), "wb") as outfile:
    pkl.dump(rxn_id_to_DFT_opted_xyz_input_ori, outfile, protocol=pkl.HIGHEST_PROTOCOL)

print("Done!")
