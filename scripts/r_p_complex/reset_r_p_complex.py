from argparse import ArgumentParser
import os
import pickle as pkl
import pandas as pd
from rdkit import Chem

from autoqm.calculation.reset_r_p_complex import reset_r_p_complex_ff_opt

parser = ArgumentParser()
parser.add_argument(
    "--input_smiles",
    type=str,
    required=True,
    help="input smiles included in a .csv file",
)
parser.add_argument(
    "--output_folder", type=str, default="output", help="output folder name"
)
parser.add_argument("--scratch_dir", type=str, required=True, help="scratch dir")
parser.add_argument(
    "--task_id",
    type=int,
    required=True,
)
parser.add_argument(
    "--num_tasks",
    type=int,
    required=True,
)
parser.add_argument(
    "--smiles_column",
    type=str,
    default="rxn_smi",
)
parser.add_argument(
    "--id_column",
    type=str,
    default="id",
)
parser.add_argument(
    "--xyz_column",
    type=str,
    default="xtb_xyz",
)

# reactant complex and product complex semiempirical optimization calculation
parser.add_argument(
    "--r_p_complex_ff_opt_folder",
    type=str,
    default="r_p_complex_ff_opt",
    help="folder for reactant complex and product complex force field optimization",
)

# specify paths
parser.add_argument(
    "--RDMC_path",
    type=str,
    required=False,
    default=None,
    help="path to RDMC to use xtb-gaussian script for xtb optimization calculation.",
)

args = parser.parse_args()

# check paths
RDMC_PATH = args.RDMC_path

assert RDMC_PATH is not None, "RDMC_PATH must be provided for ff opt"

# create directories
submit_dir = os.path.abspath(os.getcwd())
output_dir = os.path.join(submit_dir, args.output_folder)
os.makedirs(output_dir, exist_ok=True)
r_p_complex_ff_opt_dir = os.path.join(output_dir, args.r_p_complex_ff_opt_folder)
os.makedirs(r_p_complex_ff_opt_dir, exist_ok=True)
os.makedirs(args.scratch_dir, exist_ok=True)

inputs_dir = os.path.join(r_p_complex_ff_opt_dir, "inputs")
outputs_dir = os.path.join(r_p_complex_ff_opt_dir, "outputs")

# read inputs
if args.input_smiles.endswith(".csv"):
    df = pd.read_csv(args.input_smiles)
elif args.input_smiles.endswith(".pkl"):
    with open(args.input_smiles, "rb") as f:
        output = pkl.load(f)
    success_dict_list = [out[1] for out in output]
    ids = []
    rxn_smiles_list = []
    xyz_list = []

    for success_dict in success_dict_list:
        for id in success_dict:
            ids.append(id)
            rxn_smiles_list.append(success_dict[id]["rxn_smi"])
            xyz_list.append(success_dict[id]["xtb_xyz"])

    df = pd.DataFrame({"id": ids, "rxn_smi": rxn_smiles_list, "dft_xyz": xyz_list})
else:
    raise NotImplementedError
assert len(df["id"]) == len(set(df["id"])), "ids must be unique"

ts_ids = list(df[args.id_column])
rxn_smiles_list = list(df[args.smiles_column])
xyz_list = list(df[args.xyz_column])
ts_id_to_rxn_smi = dict(zip(ts_ids, rxn_smiles_list))
ts_id_to_dft_xyz = dict(zip(ts_ids, xyz_list))

print("Making inputs...")
tasks = list(zip(ts_ids, rxn_smiles_list, xyz_list))
for ts_id, rxn_smi, dft_xyz in tasks[args.task_id :: args.num_tasks]:
    ids = int(ts_id // 1000)
    suboutputs_dir = os.path.join(outputs_dir, f"outputs_{ids}")
    os.makedirs(suboutputs_dir, exist_ok=True)
    if not os.path.exists(os.path.join(suboutputs_dir, f"rxn_{ts_id}.sdf")):
        subinputs_dir = os.path.join(inputs_dir, f"inputs_{ids}")
        os.makedirs(subinputs_dir, exist_ok=True)
        if not os.path.exists(
            os.path.join(subinputs_dir, f"rxn_{ts_id}.in")
        ) and not os.path.exists(os.path.join(subinputs_dir, f"rxn_{ts_id}.tmp")):
            ts_id_input_path = os.path.join(subinputs_dir, f"rxn_{ts_id}.in")
            with open(ts_id_input_path, "w") as f:
                f.write(rxn_smi)
            print(ts_id)
            print(rxn_smi)

print("FF optimization for reactant and product complexes...")
for _ in range(5):
    for subinputs_folder in os.listdir(inputs_dir):
        ids = subinputs_folder.split("_")[1]
        subinputs_dir = os.path.join(inputs_dir, subinputs_folder)
        suboutputs_dir = os.path.join(outputs_dir, f"outputs_{ids}")
        for input_file in os.listdir(subinputs_dir):
            if ".in" in input_file:
                ts_id = int(input_file.split(".in")[0].split("rxn_")[1])
                try:
                    os.rename(
                        os.path.join(subinputs_dir, input_file),
                        os.path.join(subinputs_dir, f"rxn_{ts_id}.tmp"),
                    )
                except:
                    continue
                else:
                    rxn_smi = ts_id_to_rxn_smi[ts_id]
                    dft_xyz = ts_id_to_dft_xyz[ts_id]
                    print(ts_id)
                    print(rxn_smi)
                    reset_r_p_complex_ff_opt(
                        rxn_smi,
                        dft_xyz,
                        ts_id,
                        subinputs_dir,
                        suboutputs_dir,
                        args.scratch_dir,
                    )

print("Done!")
