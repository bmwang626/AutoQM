#!/usr/bin/env python
# coding: utf-8
import os
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
mol_ids = df["id"].tolist()
smiles_list = df["smi"].tolist()

log_paths = []
for mol_id in mol_ids:
    ids = mol_id // 1000
    log_path = os.path.join(
        "output",
        "DFT_opt_freq",
        "outputs",
        f"outputs_{ids}",
        f"{mol_id}",
        f"{mol_id}.log",
    )
    log_paths.append(log_path)

out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(
    delayed(dft_opt_freq_parser)(path, is_ts=False, check_connectivity=True, smi=smi) for path, smi in tqdm(zip(log_paths, smiles_list)) # not able to use check_connectivity=True for TS
)

failed_jobs = dict()
valid_jobs = dict()
for mol_id, (failed_job, valid_job) in zip(mol_ids, out):
    if failed_job:
        failed_jobs[mol_id] = failed_job
    if valid_job:
        valid_jobs[mol_id] = valid_job

with open(os.path.join(f"{output_file_name}.pkl"), "wb") as outfile:
    pkl.dump(valid_jobs, outfile, protocol=pkl.HIGHEST_PROTOCOL)

with open(os.path.join(f"{output_file_name}_failed.pkl"), "wb") as outfile:
    pkl.dump(failed_jobs, outfile, protocol=pkl.HIGHEST_PROTOCOL)

print(f"Total number of molecules: {len(mol_ids)}")
print(f"Total number of failed molecules: {len(failed_jobs)}")
print(failed_jobs)

mol_id_to_DFT_opted_xyz_std_ori = {}
for mol_id, valid_job in valid_jobs.items():
    mol_id_to_DFT_opted_xyz_std_ori[mol_id] = valid_job["dft_xyz_std_ori"]

with open(os.path.join(f"{output_file_name}_xyz_std_ori.pkl"), "wb") as outfile:
    pkl.dump(mol_id_to_DFT_opted_xyz_std_ori, outfile, protocol=pkl.HIGHEST_PROTOCOL)

mol_id_to_DFT_opted_xyz_input_ori = {}
for mol_id, valid_job in valid_jobs.items():
    mol_id_to_DFT_opted_xyz_input_ori[mol_id] = valid_job["dft_xyz_input_ori"]

with open(os.path.join(f"{output_file_name}_xyz_input_ori.pkl"), "wb") as outfile:
    pkl.dump(mol_id_to_DFT_opted_xyz_input_ori, outfile, protocol=pkl.HIGHEST_PROTOCOL)

print("Done!")
