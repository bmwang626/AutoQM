import os
import pandas as pd
from rdkit import Chem
from argparse import ArgumentParser

from rdmc.mol import RDKitMol

from radical_workflow.calculation.dft_calculation import dft_scf_opt

parser = ArgumentParser()
parser.add_argument(
    "--input_smiles",
    type=str,
    required=True,
    help="input smiles included in a .csv file",
)
parser.add_argument(
    "--input_geometry",
    type=str,
    required=True,
    help="input geometry included in a .sdf file",
)
parser.add_argument(
    "--output_folder", type=str, default="output", help="output folder name"
)
parser.add_argument(
    "--task_id",
    type=int,
    default=0,
    help="task id for the calculation",
)
parser.add_argument(
    "--num_tasks",
    type=int,
    default=1,
    help="number of tasks for the calculation",
)

# DFT optimization and frequency calculation
parser.add_argument(
    "--DFT_opt_freq_folder",
    type=str,
    default="DFT_opt_freq",
    help="folder for DFT optimization and frequency calculation",
)
parser.add_argument(
    "--DFT_opt_freq_theory",
    type=str,
    default="#P opt=(ts,calcall,maxcycle=32,noeig,nomicro,cartesian) scf=(xqc) iop(7/33=1) iop(2/9=2000) guess=mix wb97xd/def2svp",
    help="level of theory for the DFT calculation",
)
parser.add_argument(
    "--DFT_opt_freq_n_procs",
    type=int,
    default=16,
    help="number of process for DFT calculations",
)
parser.add_argument(
    "--DFT_opt_freq_job_ram",
    type=int,
    default=62400,  # 3900*16
    help="amount of ram (MB) allocated for each DFT calculation",
)

# specify paths
parser.add_argument(
    "--G16_path",
    type=str,
    required=False,
    default=None,
    help="path to installed Gaussian 16",
)
parser.add_argument(
    "--RDMC_path",
    type=str,
    required=False,
    default=None,
    help="path to RDMC to use xtb-gaussian script for xtb optimization calculation.",
)
parser.add_argument("--scratch_dir", type=str, required=True, help="scratch directory")

args = parser.parse_args()

G16_PATH = args.G16_path
RDMC_PATH = args.RDMC_path

submit_dir = os.path.abspath(os.getcwd())
output_dir = os.path.join(submit_dir, args.output_folder)

df = pd.read_csv(args.input_smiles)
assert len(df["id"]) == len(set(df["id"])), "ids must be unique"

assert (
    G16_PATH is not None
), f"G16_PATH must be provided for DFT optimization and frequency calculation"
assert RDMC_PATH is not None, f"RDMC_PATH must be provided to read sdf files"

# create id to smile mapping
mol_ids = df["id"].tolist()
rsmi_list = df["rsmi"].tolist()
psmi_list = df["psmi"].tolist()
smiles_list = [rsmi + ">>" + psmi for rsmi, psmi in zip(rsmi_list, psmi_list)]
mol_id_to_rxn_smi = dict(zip(mol_ids, smiles_list))
rxn_smi_to_mol_id = dict(zip(smiles_list, mol_ids))
mol_id_to_charge = dict()
mol_id_to_mult = dict()
for k, v in mol_id_to_rxn_smi.items():
    try:
        rsmi, psmi = v.split(">>")
        mol = Chem.MolFromSmiles(rsmi)
    except Exception as e:
        print(f"Cannot translate smi {v} to molecule for species {k}")

    try:
        charge = Chem.GetFormalCharge(mol)
        mol_id_to_charge[k] = charge
    except Exception as e:
        print(f"Cannot determine molecular charge for species {k} with smi {v}")

    num_radical_elec = 0
    for atom in mol.GetAtoms():
        num_radical_elec += atom.GetNumRadicalElectrons()
    mol_id_to_mult[k] = num_radical_elec + 1

os.makedirs(args.scratch_dir, exist_ok=True)
mol_ids_smis = list(zip(mol_ids, smiles_list))

# create id to xyz mapping
mols = RDKitMol.FromFile(args.input_geometry, removeHs=False, sanitize=False)
mol_id_to_xyz = dict()
for mol in mols:
    rxn_smi = mol.GetProp("_Name")
    if "_" in rxn_smi:
        rxn_smi = rxn_smi.split("_")[1]
    xyz = mol.ToXYZ()
    if rxn_smi not in rxn_smi_to_mol_id:
        print(f"Cannot find TS {rxn_smi} in the input smiles file")
        continue
    mol_id = rxn_smi_to_mol_id[rxn_smi]
    mol_id_to_xyz[mol_id] = xyz

print("DFT TS opt & freq")

for _ in range(1):
    print("Making input files for TS DFT optimization and frequency calculation")

    DFT_opt_freq_dir = os.path.join(output_dir, args.DFT_opt_freq_folder)
    os.makedirs(DFT_opt_freq_dir, exist_ok=True)
    inputs_dir = os.path.join(DFT_opt_freq_dir, "inputs")
    outputs_dir = os.path.join(DFT_opt_freq_dir, "outputs")
    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)

    for mol_id, smi in mol_ids_smis[args.task_id : len(mol_ids_smis) : args.num_tasks]:
        ids = mol_id // 1000
        input_rxns_dir = os.path.join(DFT_opt_freq_dir, "inputs", f"rxns_{ids}")
        output_rxns_dir = os.path.join(DFT_opt_freq_dir, "outputs", f"rxns_{ids}")
        os.makedirs(output_rxns_dir, exist_ok=True)
        input_rxn_dir = os.path.join(input_rxns_dir, f"rxn_{mol_id}")
        output_rxn_dir = os.path.join(output_rxns_dir, f"rxn_{mol_id}")
        os.makedirs(output_rxn_dir, exist_ok=True)

        if not os.path.exists(os.path.join(output_rxn_dir, f"{mol_id}.log")):
            if mol_id in mol_id_to_xyz:
                if not os.path.exists(
                    os.path.join(input_rxn_dir, f"{mol_id}.in")
                ) and not os.path.exists(os.path.join(input_rxn_dir, f"{mol_id}.tmp")):
                    os.makedirs(input_rxns_dir, exist_ok=True)
                    os.makedirs(input_rxn_dir, exist_ok=True)
                    with open(
                        os.path.join(input_rxn_dir, f"{mol_id}.in"),
                        "w",
                    ) as f:
                        f.write("")
                    print(mol_id)

    print("Performing DFT TS optimization and frequency calculation...")

    DFT_opt_freq_theory = args.DFT_opt_freq_theory

    for _ in range(1):
        for input_rxns_folder in os.listdir(os.path.join(DFT_opt_freq_dir, "inputs")):
            ids = int(input_rxns_folder.split("_")[1])
            input_rxns_dir = os.path.join(DFT_opt_freq_dir, "inputs", f"rxns_{ids}")
            output_rxns_dir = os.path.join(DFT_opt_freq_dir, "outputs", f"rxns_{ids}")
            for input_rxn_folder in os.listdir(input_rxns_dir):
                input_rxn_dir = os.path.join(input_rxns_dir, input_rxn_folder)
                input_file_path = os.path.join(input_rxn_dir, f"{mol_id}.in")
                if os.path.exists(input_file_path):
                    mol_id = int(input_rxn_folder.split("_")[1])
                    try:
                        print(input_file_path)
                        print(os.path.join(input_rxn_dir, f"{mol_id}.tmp"))
                        os.rename(
                            input_file_path,
                            os.path.join(input_rxn_dir, f"{mol_id}.tmp"),
                        )
                    except:
                        continue
                    else:
                        rxn_smi = mol_id_to_rxn_smi[mol_id]
                        charge = mol_id_to_charge[mol_id]
                        mult = mol_id_to_mult[mol_id]
                        xyz = mol_id_to_xyz[mol_id]

                        print(mol_id)
                        print(rxn_smi)

                        dft_scf_opt(
                            mol_id,
                            xyz,
                            G16_PATH,
                            DFT_opt_freq_theory,
                            args.DFT_opt_freq_n_procs,
                            args.DFT_opt_freq_job_ram,
                            charge,
                            mult,
                            args.scratch_dir,
                            input_rxn_dir,
                            output_rxn_dir,
                        )

    print("DFT optimization and frequency calculation done.")

print("Done!")
