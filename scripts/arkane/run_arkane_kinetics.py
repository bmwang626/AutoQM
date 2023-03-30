import os
import time
import pandas as pd
from argparse import ArgumentParser

from radical_workflow.arkane.kinetics import run_arkane_kinetics

parser = ArgumentParser()
parser.add_argument('--input_smiles', type=str, required=True,
                    help='input smiles included in a .csv file')
parser.add_argument('--task_id', type=int, default=0,
                    help='task id for the calculation',)
parser.add_argument('--num_tasks', type=int, default=1,
                    help='number of tasks for the calculation',)
parser.add_argument('--scratch_dir', type=str, required=True,
                    help='scratch directory')
parser.add_argument('--RMG_path', type=str, required=True,
                    help='path to RMG-Py')

parser.add_argument('--model_chemistry', type=str, required=True,)

args = parser.parse_args()

submit_dir = os.path.abspath(os.getcwd())

print("Making output directory")
output_dir = os.path.join(submit_dir, "output")
os.makedirs(output_dir, exist_ok=True)

print("Reading input smiles")
df = pd.read_csv(args.input_smiles, index_col=0)
assert len(df['id']) == len(set(df['id'])), "ids must be unique"

mol_ids = df['id'].tolist()
if 'rxn_smi' in df.columns:
    smiles_list = df['rxn_smi'].tolist()
else:
    raise ValueError("smiles or rxn_smi column must be present")
mol_ids_smis = list(zip(mol_ids, smiles_list))
mol_id_to_smi = dict(mol_ids_smis)
mol_id_to_row_index = {mol_id: i for i, mol_id in enumerate(mol_ids)}

print("Making helper input files")
arkane_kinetics_dir = os.path.join(output_dir, "arkane_kinetics")
os.makedirs(arkane_kinetics_dir, exist_ok=True)
inputs_dir = os.path.join(arkane_kinetics_dir, "inputs")
outputs_dir = os.path.join(arkane_kinetics_dir, "outputs")
os.makedirs(inputs_dir, exist_ok=True)
os.makedirs(outputs_dir, exist_ok=True)

print("Making dummy input files")
for mol_id, smi in mol_ids_smis[args.task_id:len(mol_ids_smis):args.num_tasks]:
    ids = str(int(int(mol_id.split("id")[1])/1000))
    subinputs_dir = os.path.join(inputs_dir, f"inputs_{ids}")
    suboutputs_dir = os.path.join(outputs_dir, f"outputs_{ids}")
    os.makedirs(suboutputs_dir, exist_ok=True)
    dummy_input_path = os.path.join(subinputs_dir, f"{mol_id}.in")
    if not os.path.exists(os.path.join(suboutputs_dir, f"{mol_id}.py")):
        os.makedirs(subinputs_dir, exist_ok=True)
        if not os.path.exists(dummy_input_path) and not os.path.exists(os.path.join(subinputs_dir, f"{mol_id}.tmp")):
            with open(dummy_input_path, "w") as f:
                f.write(mol_id)
            print(mol_id)

for subinputs_folder in os.listdir(os.path.join(arkane_kinetics_dir, "inputs")):
    ids = subinputs_folder.split("_")[1]
    subinputs_dir = os.path.join(inputs_dir, subinputs_folder)
    suboutputs_dir = os.path.join(outputs_dir, f"outputs_{ids}")
    for input_file in os.listdir(subinputs_dir):
        if ".in" in input_file:
            mol_id = input_file.split(".in")[0]
            try:
                os.rename(os.path.join(subinputs_dir, input_file), os.path.join(subinputs_dir, f"{mol_id}.tmp"))
            except:
                continue
            else:
                ids = str(int(int(mol_id.split("id")[1])/1000))
                smi = mol_id_to_smi[mol_id]
                print(mol_id)
                print(smi)

                row_index = mol_id_to_row_index[mol_id]

                start_time = time.time()
                run_arkane_kinetics(mol_id, smi, row_index, df, args.model_chemistry, subinputs_dir, suboutputs_dir, args.scratch_dir, RMG_path=args.RMG_path)
                end_time = time.time()
                print(f"Time for arkane kinetics for {mol_id} is {end_time - start_time} seconds")