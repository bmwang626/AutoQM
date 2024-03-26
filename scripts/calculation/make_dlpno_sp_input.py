from argparse import ArgumentParser
import os
import pickle as pkl
import pandas as pd
import rdkit.Chem as Chem
import subprocess

from autoqm.calculation.wft_calculation import generate_dlpno_sp_input

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
parser.add_argument(
    "--smiles_column",
    type=str,
    default="smiles",
    help="column name for smiles in the input csv file",
)
parser.add_argument(
    "--xyz_column",
    type=str,
    default="xyz",
    help="column name for xyz in the input csv file",
)
parser.add_argument(
    "--id_column",
    type=str,
    default="id",
    help="column name for id in the input csv file",
)
parser.add_argument(
    "--is_rxn",
    action="store_true",
    help="whether the input is reaction smiles",
)
parser.add_argument(
    "--xyz_has_header",
    action="store_true",
    help="whether the xyz file has header",
)

# dlpno sp
parser.add_argument(
    "--DLPNO_sp_folder",
    type=str,
    required=True,
    choices=["DLPNO_sp", "DLPNO_sp_f12"],
)
parser.add_argument(
    "--DLPNO_sp_level_of_theory",
    type=str,
    required=True,
    help="level of theory for DLPNO calculation",
)
parser.add_argument(
    "--DLPNO_sp_n_procs",
    type=int,
    default=24,
    help="number of process for DLPNO calculations",
)
parser.add_argument(
    "--DLPNO_sp_job_ram",
    type=int,
    default=4000,
    help="amount of ram (MB) per core allocated for each DLPNO calculation",
)
parser.add_argument(
    "--DLPNO_sp_cutoff_heavy_atoms",
    type=int,
    default=0,
    help="Only perform DLPNO calculation for molecules with less than this number of heavy atoms.",
)
parser.add_argument(
    "--DLPNO_sp_cutoff_total_atoms",
    type=int,
    default=0,
    help="Only perform DLPNO calculation for molecules with less than this number of total atoms (include H).",
)

# specify paths
parser.add_argument(
    "--XTB_path", type=str, required=False, default=None, help="path to installed XTB"
)
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
parser.add_argument(
    "--COSMOtherm_path",
    type=str,
    required=False,
    default=None,
    help="path to COSMOthermo",
)
parser.add_argument(
    "--COSMO_database_path",
    type=str,
    required=False,
    default=None,
    help="path to COSMO_database",
)
parser.add_argument(
    "--ORCA_path", type=str, required=False, default=None, help="path to ORCA"
)

args = parser.parse_args()

XTB_PATH = args.XTB_path
G16_PATH = args.G16_path
RDMC_PATH = args.RDMC_path
COSMOTHERM_PATH = args.COSMOtherm_path
COSMO_DATABASE_PATH = args.COSMO_database_path
ORCA_PATH = args.ORCA_path
assert ORCA_PATH is not None, "ORCA_PATH must be provided for dlpno sp calc"

submit_dir = os.path.abspath(os.getcwd())
output_dir = os.path.join(submit_dir, args.output_folder)
DLPNO_sp_dir = os.path.join(output_dir, args.DLPNO_sp_folder)
os.makedirs(DLPNO_sp_dir, exist_ok=True)

# input files
df = pd.read_csv(args.input_smiles)
assert len(df[args.id_column]) == len(set(df[args.id_column])), "ids must be unique"

if args.xyz_has_header:
    xyzs = df[args.xyz_column].apply(lambda x: "\n".join(x.split("\n")[2:]))
else:
    xyzs = df[args.xyz_column]
xyz_DFT_opt_dict = dict(zip(df[args.id_column], xyzs))

mol_ids = list(df[args.id_column])
if not args.is_rxn:
    smiles_list = list(df[args.smiles_column])
else:
    smiles_list = list(df[args.smiles_column])
    smiles_list = [smi.split(">>")[0] for smi in smiles_list]  # use the reactant smiles

# create id to property mapping
mol_id_to_smi_dict = dict(zip(mol_ids, smiles_list))
mol_id_to_charge_dict = dict()
mol_id_to_mult_dict = dict()
mol_id_to_num_heavy_atoms_dict = dict()
mol_id_to_num_total_atoms_dict = dict()

# smi parsing settings
params = Chem.SmilesParserParams()
# rdkit by default removes explicit Hs, to make sure total atom count is correct we need to include all Hs in smiles
params.removeHs = False

for k, v in mol_id_to_smi_dict.items():
    try:
        mol = Chem.MolFromSmiles(v, params)
    except Exception as e:
        print(f"Cannot translate smi {v} to molecule for species {k}")

    try:
        charge = Chem.GetFormalCharge(mol)
        mol_id_to_charge_dict[k] = charge
    except Exception as e:
        print(f"Cannot determine molecular charge for species {k} with smi {v}")

    num_radical_elec = 0
    for atom in mol.GetAtoms():
        num_radical_elec += atom.GetNumRadicalElectrons()
    mol_id_to_mult_dict[k] = num_radical_elec + 1

    num_total_atoms = 0
    num_heavy_atoms = 0
    for atom in mol.GetAtoms():
        num_total_atoms += 1
        if atom.GetAtomicNum() > 1:
            num_heavy_atoms += 1
    mol_id_to_num_heavy_atoms_dict[k] = num_heavy_atoms
    mol_id_to_num_total_atoms_dict[k] = num_total_atoms

inputs_dir = os.path.join(DLPNO_sp_dir, "inputs")
os.makedirs(inputs_dir, exist_ok=True)
outputs_dir = os.path.join(DLPNO_sp_dir, "outputs")
os.makedirs(outputs_dir, exist_ok=True)

print("Make dlpno input files...")


def get_maxcore(log_path):
    proc = subprocess.run(
        f"head -n 200 {log_path} | grep maxcore", shell=True, capture_output=True
    )
    strings = proc.stdout.split()
    if strings:
        return int(strings[-1])
    return None


def tail(f, n):
    proc = subprocess.Popen(["tail", "-n", str(n), f], stdout=subprocess.PIPE)
    lines = proc.stdout.readlines()
    return lines


def has_mdci_error(log_path, bypass=False, exclude_convergence=True):
    """
    Check MDCI error in ORCA log. Most MDCI errors are due to ram issue or unexpected kill of the job. Use bypass to skip this test when needed as this might interference with maxcore check.
    """

    if bypass:
        return False

    lines = tail(log_path, 30)

    # convergence error might also shown as MDCI error
    if exclude_convergence:
        for line in lines:
            if b"ORCA finished by error termination in SCF" in line:
                return False
            if b"This wavefunction IS NOT CONVERGED!" in line:
                return False

    for line in lines:
        if b"ORCA finished by error termination in MDCI" in line:
            return True
    else:
        return False


def has_max_core_error(log_path, current_maxcore):
    lines = tail(log_path, 30)
    for line in lines:
        if (
            (b"Please increase MaxCore - Skipping calculation" in line)
            or (b"ORCA finished by error termination in GTOInt" in line)
            or (b"ORCA finished by error termination in MDCI" in line)
        ):
            maxcore = get_maxcore(log_path)
            if maxcore is not None and maxcore < current_maxcore:
                return True
    # max core error can show up as SCF error
    scf_error = False
    wave_function_error = False
    for line in lines:
        if b"ORCA finished by error termination in SCF" in line:
            scf_error = True
        if b"This wavefunction IS NOT CONVERGED!" in line:
            wave_function_error = True
    if scf_error and not wave_function_error:
        return True

    return False


def has_wave_function_error(log_path):
    lines = tail(log_path, 30)
    # find wave function error
    wave_function_error = False
    for line in lines:
        if b"This wavefunction IS NOT CONVERGED!" in line:
            wave_function_error = True
    if wave_function_error:
        return True

    return False


# added by oscar Jan 3, 2023
def has_coordinates_error(log_path):
    lines = tail(log_path, 30)
    # find coordinates error
    cord_error = False
    for line in lines:
        if b"You must have a" in line:
            cord_error = True
        if b"block in your input" in line:
            cord_error = True
    if cord_error:
        return True

    return False


mol_ids_smis = list(zip(mol_ids, smiles_list))
for mol_id, smi in mol_ids_smis[args.task_id :: args.num_tasks]:

    # determine if we need to skip some jobs based on atom counts
    if not any([args.DLPNO_sp_cutoff_heavy_atoms, args.DLPNO_sp_cutoff_total_atoms]):
        # both keywords are not specified, proceed to calculate all
        pass
    elif all([args.DLPNO_sp_cutoff_heavy_atoms, args.DLPNO_sp_cutoff_total_atoms]):
        # both keywords are specified, check all conditions
        if mol_id_to_num_heavy_atoms_dict[mol_id] > args.DLPNO_sp_cutoff_heavy_atoms:
            continue
        if mol_id_to_num_total_atoms_dict[mol_id] > args.DLPNO_sp_cutoff_total_atoms:
            continue
        # cases for only one keyword is specified
    elif args.DLPNO_sp_cutoff_heavy_atoms:
        # skip based on heavy atom count
        if mol_id_to_num_heavy_atoms_dict[mol_id] > args.DLPNO_sp_cutoff_heavy_atoms:
            continue
    elif args.DLPNO_sp_cutoff_total_atoms:
        if mol_id_to_num_total_atoms_dict[mol_id] > args.DLPNO_sp_cutoff_total_atoms:
            continue

    ids = mol_id // 1000
    subinputs_dir = os.path.join(inputs_dir, f"inputs_{ids}")
    suboutputs_dir = os.path.join(outputs_dir, f"outputs_{ids}")
    os.makedirs(suboutputs_dir, exist_ok=True)
    log_path = os.path.join(suboutputs_dir, f"{mol_id}.log")
    DLPNO_sp_level_of_theory = args.DLPNO_sp_level_of_theory
    if mol_id in xyz_DFT_opt_dict:
        if os.path.exists(log_path):
            # check if maxcore error
            # if has_max_core_error(log_path, args.DLPNO_sp_job_ram):
            #     print(f"maxcore error for {mol_id}, removing...")
            #     try:
            #         os.remove(log_path)
            #     except FileNotFoundError:
            #         print(f"file {log_path} not found, already removed?")

            # check if mdci error, set bypass if want to skip (not every restart need this, useful for core/ram change mostly)
            if has_mdci_error(log_path, bypass=False):
                print(f"mdci error for {mol_id}, removing...")
                try:
                    os.remove(log_path)
                except FileNotFoundError:
                    print(f"file {log_path} not found, already removed?")

            # check for cord error
            if has_coordinates_error(log_path):
                print(f"cords error for {mol_id}, removing...")
                try:
                    os.remove(log_path)
                except FileNotFoundError:
                    print(f"file {log_path} not found, already removed?")

            # check if wave function error
        #     if has_wave_function_error(log_path):
        #         try:
        #             if args.DLPNO_sp_level_of_theory == "uHF UNO DLPNO-CCSD(T)-F12D cc-pvtz-f12 def2/J cc-pvqz/c cc-pvqz-f12-cabs RIJCOSX VeryTightSCF NormalPNO":
        #                 DLPNO_sp_level_of_theory = "uHF UNO DLPNO-CCSD(T)-F12D cc-pvtz-f12 def2/J cc-pvqz/c cc-pvqz-f12-cabs RIJCOSX NormalSCF NormalPNO"
        #                 print(f"wave function error for {mol_id}, removing...")
        #                 os.remove(log_path)
        #             else:
        #                 print(f"wave function error for {mol_id}, already using NormalSCF, skipping...")
        #         except FileNotFoundError:
        #             print(f"file {log_path} not found, already removed?")
        if not os.path.exists(log_path):
            os.makedirs(subinputs_dir, exist_ok=True)
            mol_id_path = os.path.join(subinputs_dir, f"{mol_id}.in")
            if not os.path.exists(
                os.path.join(subinputs_dir, f"{mol_id}.tmp")
            ) and not os.path.exists(mol_id_path):
                charge = mol_id_to_charge_dict[mol_id]
                mult = mol_id_to_mult_dict[mol_id]
                coords = xyz_DFT_opt_dict[mol_id].strip()
                script = generate_dlpno_sp_input(
                    DLPNO_sp_level_of_theory,
                    coords,
                    charge,
                    mult,
                    args.DLPNO_sp_job_ram,
                    args.DLPNO_sp_n_procs,
                )

                print(f"Generating input file for {mol_id}...")
                with open(mol_id_path, "w+") as f:
                    f.write(script)
    else:
        print(f"Cannot find xyz for {mol_id}")
        if os.path.exists(log_path):
            print(f"Removing {log_path}...")
            os.remove(log_path)

print("Done!")
