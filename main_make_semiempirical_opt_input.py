from argparse import ArgumentParser
import os
import pickle as pkl
import pandas as pd

parser = ArgumentParser()
parser.add_argument('--input_smiles', type=str, required=True,
                    help='input smiles included in a .csv file')
parser.add_argument('--output_folder', type=str, default='output',
                    help='output folder name')
parser.add_argument('--xyz_FF_dict', type=str, required=True,
                    help='pickled dict mapping from mol_id to confs xyz')
parser.add_argument('--task_id', type=int, default=0,
                    help='task id for the calculation',)
parser.add_argument('--num_tasks', type=int, default=1,
                    help='number of tasks for the calculation',)

# semiempirical optimization calculation
parser.add_argument('--semiempirical_opt_folder', type=str, default='semiempirical_opt',
                    help='folder for semiempirical optimization')
parser.add_argument('--semiempirical_method', type=str, default='all',
                    help='method used for semiempirical optimization. Options are GFN2-XTB, pm7, am1, and all. If all is chosen, GFN2-XTB will be used first, pm7 second, and am1 third if the previous method does not work.')
parser.add_argument('--gaussian_semiempirical_opt_theory', type=str, default='#opt=(calcall,maxcycle=128,noeig,nomicro,cartesian)',
                    help='level of theory for the Gaussian semiempirical calculation')
parser.add_argument('--gaussian_semiempirical_opt_n_procs', type=int, default=4,
                    help='number of process for Gaussian semiempirical calculations')
parser.add_argument('--gaussian_semiempirical_opt_job_ram', type=int, default=3000,
                    help='amount of ram (MB) allocated for each Gaussian semiempirical calculation')

# specify paths
parser.add_argument('--XTB_path', type=str, required=False, default=None,
                    help='path to installed XTB')
parser.add_argument('--G16_path', type=str, required=False, default=None,
                    help='path to installed Gaussian 16')
parser.add_argument('--RDMC_path', type=str, required=False, default=None,
                    help='path to RDMC to use xtb-gaussian script for xtb optimization calculation.')
parser.add_argument('--COSMOtherm_path', type=str, required=False, default=None,
                    help='path to COSMOthermo')
parser.add_argument('--COSMO_database_path', type=str, required=False, default=None,
                    help='path to COSMO_database')
parser.add_argument('--ORCA_path', type=str, required=False, default=None,
                    help='path to ORCA')

args = parser.parse_args()

XTB_PATH = args.XTB_path
G16_PATH = args.G16_path
RDMC_PATH = args.RDMC_path
COSMOTHERM_PATH = args.COSMOtherm_path
COSMO_DATABASE_PATH = args.COSMO_database_path
ORCA_PATH = args.ORCA_path

submit_dir = os.path.abspath(os.getcwd())
output_dir = os.path.join(submit_dir, args.output_folder)
semiempirical_opt_dir = os.path.join(output_dir, args.semiempirical_opt_folder)
os.makedirs(semiempirical_opt_dir, exist_ok=True)

df = pd.read_csv(args.input_smiles, index_col=0)
assert len(df['id']) == len(set(df['id'])), "ids must be unique"

with open(args.xyz_FF_dict, "rb") as f:
    xyz_FF_dict = pkl.load(f)

assert XTB_PATH is not None, "XTB_PATH must be provided for semiempirical opt"
assert G16_PATH is not None, "G16_PATH must be provided for semiempirical opt"
assert RDMC_PATH is not None, "RDMC_PATH must be provided for semiempirical opt"

mol_ids = list(df["id"])
smiles_list = list(df["smiles"])
inputs_dir = os.path.join(semiempirical_opt_dir, "inputs")
os.makedirs(inputs_dir, exist_ok=True)
outputs_dir = os.path.join(semiempirical_opt_dir, "outputs")
os.makedirs(outputs_dir, exist_ok=True)

mol_ids_smis = zip(mol_ids, smiles_list)
for mol_id, smi in mol_ids_smis[args.task_id:args.num_tasks:len(mol_ids_smis)]:
    if mol_id in xyz_FF_dict:
        ids = str(int(int(mol_id.split("id")[1])/1000))
        subinputs_dir = os.path.join(inputs_dir, f"inputs_{ids}")
        suboutputs_dir = os.path.join(outputs_dir, f"outputs_{ids}")
        os.makedirs(suboutputs_dir, exist_ok=True)
        mol_id_path = os.path.join(subinputs_dir, f"{mol_id}.in")
        if not os.path.exists(os.path.join(suboutputs_dir, f"{mol_id}.tar")) and not os.path.exists(mol_id_path):
            os.makedirs(subinputs_dir, exist_ok=True)
            with open(mol_id_path, "w") as f:
                f.write(mol_id)
            print(mol_id)
