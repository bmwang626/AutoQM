import pickle as pkl
import pandas as pd
import sys

# Load input
input_path = sys.argv[1]
output_name = sys.argv[2]

input_df = pd.read_csv(input_path)
mol_id_to_index = {mol_id: i for i, mol_id in enumerate(input_df["id"])}

# check FF results
with open(f"{output_name}_ff_opted_results_failed.pkl", "rb") as f:
    failed_jobs = pkl.load(f)

failed_reason = [""] * len(input_df)
for mol_id in failed_jobs:
    index = mol_id_to_index[mol_id]
    failed_reason[index] = failed_jobs[mol_id].get("reason", "")

input_df["ff_failed_reason"] = failed_reason

# check semiempirical results
with open(f"{output_name}_semiempirical_opted_results_failed.pkl", "rb") as f:
    failed_jobs = pkl.load(f)

failed_reason = [""] * len(input_df)
for mol_id in failed_jobs:
    index = mol_id_to_index[mol_id]
    failed_reason[index] = failed_jobs[mol_id].get("reason", "")

input_df["semiempirical_failed_reason"] = failed_reason

# check DFT results
with open(f"{output_name}_dft_opted_results_failed.pkl", "rb") as f:
    failed_jobs = pkl.load(f)

failed_reason = [""] * len(input_df)
for mol_id in failed_jobs:
    index = mol_id_to_index[mol_id]
    failed_reason[index] = failed_jobs[mol_id].get("reason", "")

input_df["dft_failed_reason"] = failed_reason

# check dlpno results
with open(f"{output_name}_dlpno_sp_results_failed.pkl", "rb") as f:
    failed_jobs = pkl.load(f)

failed_reason = [""] * len(input_df)
for mol_id in failed_jobs:
    index = mol_id_to_index[mol_id]
    failed_reason[index] = failed_jobs[mol_id].get("reason", "")

input_df["dlpno_failed_reason"] = failed_reason

# check COSMO results
with open(f"{output_name}_cosmo_results_failed.pkl", "rb") as f:
    failed_jobs = pkl.load(f)

failed_reason = [""] * len(input_df)
for mol_id in failed_jobs:
    index = mol_id_to_index[mol_id]
    failed_reason[index] = "file not found"

input_df["cosmo_failed_reason"] = failed_reason

input_df.to_csv(f"{output_name}_failed_reasons.csv", index=False)
