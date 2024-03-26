import os
import pickle as pkl
import shutil
import tarfile

import pandas as pd
from autoqm.calculation.utils import REPLACE_LETTER
from rdkit import Chem

REPLACE_LETTER = {"(": "_", ")": "_", "'": "_"}

def mol2xyz(mol):
    return mol.ToXYZ()

def mol2charge(mol):
    return mol.GetFormalCharge()

def mol2mult(mol):
    num_radical_elec = 0
    for atom in mol.GetAtoms():
        num_radical_elec += atom.GetNumRadicalElectrons()
    return num_radical_elec + 1

def add_shared_arguments(parser):
    input_parser = parser.add_argument_group('Input')
    input_parser.add_argument('--input_file', type=str, required=True,
                        help='input CSV file containing the species information')
    input_parser.add_argument('--id_column', type=str, default='id',
                        help='column name for the species id')
    input_parser.add_argument('--smiles_column', type=str, default='smiles',
                        help='column name for the SMILES string')
    input_parser.add_argument('--xyz_column', type=str, default=None,
                        help='column name for the xyz string')
    
    parser.add_argument('--scratch_dir', type=str, required=True,
                        help='scfratch directory')
    parser.add_argument('--xyz_DFT_opt_dict', type=str, default=None,
                        help='pickle file containing a dictionary to map between the mol_id and DFT-optimized xyz for following calculations',)
    parser.add_argument('--task_id', type=int, default=0,
                        help='task id for the calculation',)
    parser.add_argument('--num_tasks', type=int, default=1,
                        help='number of tasks for the calculation',)
    parser.add_argument('--RDMC_path', type=str, required=False, default=None,
                        help='path to RDMC to use xtb-gaussian script for xtb optimization calculation.')
    return parser

def add_cosmo_arguments(parser):
    # Turbomole and COSMO calculation
    parser.add_argument('--COSMO_folder', type=str, default='COSMO_calc',
                        help='folder for COSMO calculation',)
    parser.add_argument('--COSMO_temperatures', type=str, nargs="+", required=False, default=['297.15', '298.15', '299.15'],
                        help='temperatures used for COSMO calculation')
    parser.add_argument('--COSMO_input_pure_solvents', type=str, required=False, default='common_solvent_list_final.csv',
                        help='input file containing pure solvents used for COSMO calculation.')
    parser.add_argument('--COSMOtherm_path', type=str, required=False, default=None,
                        help='path to COSMOthermo')
    parser.add_argument('--COSMO_database_path', type=str, required=False, default=None,
                        help='path to COSMO_database')

def add_xtb_arguments(parser):
    parser.add_argument('--XTB_path', type=str, required=False, default=None,
                        help='path to installed XTB')
    parser.add_argument('--G16_path', type=str, required=False, default=None,
                        help='path to installed Gaussian 16')
    return parser

def add_dlpno_arguments(parser):
    parser.add_argument('--ORCA_path', type=str, required=False, default=None,
                        help='path to ORCA')
    return parser

def add_qm_des_arguments(parser):
    qm_des_parser = parser.add_argument_group('QM descriptor calculation')
    qm_des_parser.add_argument('--title_card', type=str, required=True,
                    help='level of theory for QM descriptor calculation')
    qm_des_parser.add_argument('--G16_path', type=str, required=True,
                    help='path to installed Gaussian 16')
    return parser



df = pd.read_csv(args.input_smiles)
mol_ids = list(df.id)

if args.xyz_DFT_opt_dict:
    # input files
    with open(args.xyz_DFT_opt_dict, "rb") as f:
        xyz_DFT_opt_dict = pkl.load(f)
else:
    xyz_DFT_opt_dict = dict(zip(df.id, df.std_xyz_str.str.split("\n\n").str[1]))

if "smiles" in df.columns:
    mol_smis = list(df.smiles)
elif "rxn_smi" in df.columns:
    mol_smis = list(df.rxn_smi.str.split(">>").str[0])
elif "smi" in df.columns:
    mol_smis = list(df.smi)
else:
    raise ValueError("Cannot find smiles or rxn_smi in input file.")

mol_id_to_smi_dict = dict(zip(mol_ids, mol_smis))
mol_id_to_charge_dict = dict()
mol_id_to_mult_dict = dict()
for k, v in mol_id_to_smi_dict.items():
    try:
        mol = Chem.MolFromSmiles(v)
    except Exception as e:
        print(f'Cannot translate smi {v} to molecule for species {k}')

    try:
        charge = Chem.GetFormalCharge(mol)
        mol_id_to_charge_dict[k] = charge
    except Exception as e:
        print(f'Cannot determine molecular charge for species {k} with smi {v}')

    num_radical_elec = 0
    for atom in mol.GetAtoms():
        num_radical_elec += atom.GetNumRadicalElectrons()
    mol_id_to_mult_dict[k] =  num_radical_elec + 1

submit_dir = os.path.abspath(os.getcwd())
project_dir = os.path.abspath(os.path.join(args.output_folder))
COSMO_dir = os.path.join(project_dir, args.COSMO_folder)

df_pure = pd.read_csv(os.path.join(submit_dir,args.COSMO_input_pure_solvents))
df_pure = df_pure.reset_index()
last_cosmo_name = df_pure.loc[len(df_pure.index)-1, "cosmo_name"]
last_cosmo_name_replaced = "".join(letter if letter not in REPLACE_LETTER else REPLACE_LETTER[letter] for letter in last_cosmo_name)
COSMOTHERM_PATH = args.COSMOtherm_path
COSMO_DATABASE_PATH = args.COSMO_database_path
assert COSMOTHERM_PATH is not None and COSMO_DATABASE_PATH is not None, "COSMOTHERM_PATH and COSMO_DATABASE_PATH must be provided for COSMO calc"

print("Making inputs and outputs dir...")
inputs_dir = os.path.join(COSMO_dir, "inputs")
os.makedirs(inputs_dir, exist_ok=True)
outputs_dir = os.path.join(COSMO_dir, "outputs")
os.makedirs(outputs_dir, exist_ok=True)

print("Making helper input files...")
print(f"Task id: {args.task_id}")
print(f"Number of tasks: {args.num_tasks}")

mol_ids_smis = list(zip(mol_ids, mol_smis))

for mol_id, smi in mol_ids_smis[args.task_id::args.num_tasks]:
    if mol_id in xyz_DFT_opt_dict:
        print(f"Mol id: {mol_id} in xyz dict")

        ids = mol_id // 1000
        subinputs_dir = os.path.join(inputs_dir, f"inputs_{ids}")
        suboutputs_dir = os.path.join(outputs_dir, f"outputs_{ids}")
        os.makedirs(suboutputs_dir, exist_ok=True)
        input_file_path = os.path.join(subinputs_dir, f"{mol_id}.in")
        tmp_input_file_path = os.path.join(subinputs_dir, f"{mol_id}.tmp")
        tar_file_path = os.path.join(suboutputs_dir, f"{mol_id}.tar")
        mol_tmp_dir = os.path.join(suboutputs_dir, f"{mol_id}")
        mol_tmp_log_path = os.path.join(mol_tmp_dir, f"{mol_id}.log")

        if os.path.exists(tar_file_path):
            tar = tarfile.open(tar_file_path, "r")
            member_basename_list = set(os.path.basename(member.name) for member in tar.getmembers())
            if any(f"_{last_cosmo_name_replaced}.tab" in member_basename for member_basename in member_basename_list):
                tar.close()
                print(f"COSMO-RS calculation for {mol_id} already finished.")

                if os.path.exists(mol_tmp_dir):
                    shutil.rmtree(mol_tmp_dir)

                continue
            tar.close()

        if os.path.exists(mol_tmp_log_path):

            with open(mol_tmp_log_path, "r") as f:
                lines = f.readlines()
                print(lines)
                if any("calculation did not converge" in line for line in lines):
                    print(mol_tmp_log_path)
                    print(f"Turbomole calculation for {mol_id} did not converge. Skipping...")

                    continue

            
        os.makedirs(subinputs_dir, exist_ok=True)
        if not os.path.exists(input_file_path) and not os.path.exists(tmp_input_file_path):
            with open(input_file_path, "w+") as f:
                f.write("")
            print(f"Making helper input file for {mol_id}...")

    else:
        print(f"Mol id: {mol_id} not in xyz dict")

print("Starting COSMO calculations...")
for subinputs_folder in os.listdir(os.path.join(COSMO_dir, "inputs")):
    ids = int(subinputs_folder.split("_")[1])
    subinputs_dir = os.path.join(COSMO_dir, "inputs", subinputs_folder)
    suboutputs_dir = os.path.join(COSMO_dir, "outputs", f"outputs_{ids}")
    for input_file in os.listdir(subinputs_dir):
        if input_file.endswith(".in"):
            input_file_path = os.path.join(subinputs_dir, input_file)
            mol_id = int(input_file.split(".in")[0])
            tmp_input_file_path = os.path.join(subinputs_dir, f"{mol_id}.tmp")
            if not os.path.exists(tmp_input_file_path):
                try:
                    os.rename(input_file_path, tmp_input_file_path)
                except:
                    continue
                else:
                    print(f"Starting COSMO-RS and Turbomole calculation for {mol_id}...")
                    ids = mol_id // 1000
                    charge = mol_id_to_charge_dict[mol_id]
                    mult = mol_id_to_mult_dict[mol_id]
                    coords = xyz_DFT_opt_dict[mol_id]
                    tmp_mol_dir = os.path.join(suboutputs_dir, f"{mol_id}")
                    os.makedirs(tmp_mol_dir, exist_ok=True)
                    cosmo_calc(mol_id, COSMOTHERM_PATH, COSMO_DATABASE_PATH, charge, mult, args.COSMO_temperatures, df_pure, coords, args.scratch_dir, tmp_mol_dir, suboutputs_dir, subinputs_dir)
                    print(f"Finished COSMO-RS and Turbomole calculation for {mol_id}")

print("Done!")