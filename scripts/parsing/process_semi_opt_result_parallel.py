#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import pickle as pkl
from tqdm import tqdm
from joblib import Parallel, delayed

from autoqm.parser.semiempirical_opt_parser import semiempirical_opt_parser

input_smiles_path = sys.argv[1]
output_file_name = sys.argv[2]
n_jobs = int(sys.argv[3])

##
# input_smiles_path = "inputs/reactants_products_aug11b_inputs.csv"
# output_file_name = "reactants_products_aug11b"
# n_jobs = 1

df = pd.read_csv(input_smiles_path)
mol_id_to_smi = dict(zip(df.id, df.smiles))
mol_ids = list(df.id)

##
# mol_ids = mol_ids[:500]

out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(delayed(semiempirical_opt_parser)(mol_id, mol_id_to_smi[mol_id]) for mol_id in tqdm(mol_ids))

failed_jobs = dict()
valid_jobs = dict()
for failed_job, valid_job in out:
    failed_jobs.update(failed_job)
    valid_jobs.update(valid_job)

with open(os.path.join(f'{output_file_name}.pkl'), 'wb') as outfile:
    pkl.dump(valid_jobs, outfile, protocol=pkl.HIGHEST_PROTOCOL)

with open(os.path.join(f'{output_file_name}_failed.pkl'), 'wb') as outfile:
    pkl.dump(failed_jobs, outfile, protocol=pkl.HIGHEST_PROTOCOL)

print(f"Total number of molecules: {len(mol_ids)}")
all_failed_jobs = [job for job in failed_jobs.values() if "reason" in job]
print(f"Total number of all failed molecules: {len(all_failed_jobs)}")
print(all_failed_jobs)
