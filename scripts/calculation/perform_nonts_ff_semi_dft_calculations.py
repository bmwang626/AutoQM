from argparse import ArgumentParser
import os
import pandas as pd
import time

from rdkit import Chem
from rdmc.mol import RDKitMol

from autoqm.calculation.ff_conf_generation import _genConf
from autoqm.calculation.semiempirical_calculation import semiempirical_opt
from autoqm.calculation.dft_calculation import dft_scf_opt
from autoqm.parser.semiempirical_opt_parser import semiempirical_opt_parser, get_mol_id_to_semiempirical_opted_xyz

parser = ArgumentParser()
parser.add_argument('--input_smiles', type=str, required=True,
                    help='input smiles included in a .csv file')
parser.add_argument('--output_folder', type=str, default='output',
                    help='output folder name')
parser.add_argument('--task_id', type=int, default=0,
                    help='task id for the calculation',)
parser.add_argument('--num_tasks', type=int, default=1,
                    help='number of tasks for the calculation',)

# conformer searching
parser.add_argument('--FF_conf_folder', type=str, default='FF_conf',
                    help='Folder name for FF searched conformers')
parser.add_argument('--max_n_conf', type=int, default=800,
                    help='maximum number of FF conformers. nc = 3**n_rotatable_bonds, n_conf = nc if nc < max_n_conf else max_n_conf')
parser.add_argument('-max_conf_try', type=int, default=2000,
                    help='maximum attempt for conformer generating, '
                         'this is useful for molecules with many chiral centers.')
parser.add_argument('-rmspre', type=float, required=False, default=0.1,
                        help='rms threshold pre optimization')
parser.add_argument('--rmspost', type=float, required=False, default=0.4,
                    help='rms threshold post FF minimization')
parser.add_argument('--E_cutoff_fraction', type=float, required=False, default=0.5,
                    help='energy window for FF minimization.')
parser.add_argument('--n_lowest_E_confs_to_save', type=int, default=10,
                    help='number of lowest energy conformers to save')

# semiempirical optimization calculation
parser.add_argument('--semiempirical_opt_folder', type=str, default='semiempirical_opt',
                    help='folder for semiempirical optimization')
parser.add_argument('--gaussian_semiempirical_opt_theory', type=str, default='#opt=(calcall,maxcycle=128,noeig,nomicro,cartesian)',
                    help='level of theory for the Gaussian semiempirical calculation')
parser.add_argument('--gaussian_semiempirical_opt_n_procs', type=int, default=16,
                    help='number of process for Gaussian semiempirical calculations')
parser.add_argument('--gaussian_semiempirical_opt_job_ram', type=int, default=8000,
                    help='amount of ram (MB) allocated for each Gaussian semiempirical calculation')

# DFT optimization and frequency calculation
parser.add_argument('--DFT_opt_freq_folder', type=str, default='DFT_opt_freq',
                    help='folder for DFT optimization and frequency calculation',)
parser.add_argument('--DFT_opt_freq_theory', type=str, default='#P opt=(calcfc,maxcycle=128,noeig,nomicro,cartesian) freq scf=(xqc) iop(7/33=1) iop(2/9=2000) guess=mix wb97xd/def2svp',
                    help='level of theory for the DFT calculation')
parser.add_argument('--DFT_opt_freq_theory_backup', type=str, default='#P opt=(calcall,maxcycle=64,noeig,nomicro,cartesian) freq scf=(tight, xqc) iop(7/33=1) iop(2/9=2000) guess=mix wb97xd/def2svp',
                    help='level of theory for the DFT calculation if DFT_opt_freq_theory failed')
parser.add_argument('--DFT_opt_freq_n_procs', type=int, default=16,
                    help='number of process for DFT calculations')
parser.add_argument('--DFT_opt_freq_job_ram', type=int, default=62400, #3900*16
                    help='amount of ram (MB) allocated for each DFT calculation')

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
parser.add_argument('--scratch_dir', type=str, required=True,
                    help='scratch directory')

args = parser.parse_args()

XTB_PATH = args.XTB_path
G16_PATH = args.G16_path
RDMC_PATH = args.RDMC_path
COSMOTHERM_PATH = args.COSMOtherm_path
COSMO_DATABASE_PATH = args.COSMO_database_path
ORCA_PATH = args.ORCA_path

submit_dir = os.path.abspath(os.getcwd())
output_dir = os.path.join(submit_dir, args.output_folder)

df = pd.read_csv(args.input_smiles)
assert len(df['id']) == len(set(df['id'])), "ids must be unique"

assert XTB_PATH is not None, f"XTB_PATH must be provided for GFNFF conformer search"
assert G16_PATH is not None, f"G16_PATH must be provided for semiempirical optimization and DFT optimization and frequency calculation"
assert RDMC_PATH is not None, f"RDMC_PATH must be provided for xtb optimization calculation"

# create id to smile mapping
mol_ids = df['id'].tolist()
smiles_list = df['smi'].tolist()
mol_id_to_smi = dict(zip(mol_ids, smiles_list))
mol_id_to_charge = dict()
mol_id_to_mult = dict()
for k, v in mol_id_to_smi.items():
    try:
        mol = Chem.MolFromSmiles(v)
    except Exception as e:
        print(f'Cannot translate smi {v} to molecule for species {k}')

    try:
        charge = Chem.GetFormalCharge(mol)
        mol_id_to_charge[k] = charge
    except Exception as e:
        print(f'Cannot determine molecular charge for species {k} with smi {v}')

    num_radical_elec = 0
    for atom in mol.GetAtoms():
        num_radical_elec += atom.GetNumRadicalElectrons()
    mol_id_to_mult[k] =  num_radical_elec + 1

os.makedirs(args.scratch_dir, exist_ok=True)
mol_ids_smis = list(zip(mol_ids, smiles_list))

print("Force-field conformer search -> semiempirical optimization -> DFT optimization & frequency calculation")

for _ in range(1):

    print("Making input files for conformer searching...")

    FF_conf_dir = os.path.join(output_dir, args.FF_conf_folder)
    os.makedirs(FF_conf_dir, exist_ok=True)
    inputs_dir = os.path.join(FF_conf_dir, "inputs")
    outputs_dir = os.path.join(FF_conf_dir, "outputs")
    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)

    for mol_id, smi in mol_ids_smis[args.task_id:len(mol_ids_smis):args.num_tasks]:
        ids = mol_id // 1000
        subinputs_dir = os.path.join(FF_conf_dir, "inputs", f"inputs_{ids}")
        suboutputs_dir = os.path.join(FF_conf_dir, "outputs", f"outputs_{ids}")
        os.makedirs(suboutputs_dir, exist_ok=True)
        mol_id_path = os.path.join(subinputs_dir, f"{mol_id}.in")
        if not os.path.exists(os.path.join(suboutputs_dir, f"{mol_id}_confs.sdf")):
            os.makedirs(subinputs_dir, exist_ok=True)
            if not os.path.exists(mol_id_path) and not os.path.exists(os.path.join(subinputs_dir, f"{mol_id}.tmp")):
                with open(mol_id_path, "w") as f:
                    f.write("")
                print(f"Making input file for conformer searching for {mol_id}...")

    print("Conformer searching with force field...")

    conf_search_FFs = ["GFNFF", "MMFF94s"]
    for _ in range(1):
        for subinputs_folder in os.listdir(os.path.join(FF_conf_dir, "inputs")):
            ids = subinputs_folder.split("_")[1]
            subinputs_dir = os.path.join(FF_conf_dir, "inputs", subinputs_folder)
            suboutputs_dir = os.path.join(FF_conf_dir, "outputs", f"outputs_{ids}")
            for input_file in os.listdir(subinputs_dir):
                if ".in" in input_file:
                    mol_id = int(input_file.split(".in")[0])
                    try:
                        os.rename(os.path.join(subinputs_dir, input_file), os.path.join(subinputs_dir, f"{mol_id}.tmp"))
                    except:
                        continue
                    else:
                        ids = mol_id // 1000
                        smi = mol_id_to_smi[mol_id]
                        print(f"Conformer searching with force field for {mol_id} {smi}...")
                        start_time = time.time()
                        _genConf(smi, mol_id, XTB_PATH, conf_search_FFs, args.max_n_conf, args.max_conf_try, args.rmspre, args.E_cutoff_fraction, args.rmspost, args.n_lowest_E_confs_to_save, args.scratch_dir, suboutputs_dir, subinputs_dir)
                        end_time = time.time()
                        print(f"Time for conformer search for {mol_id} took {end_time - start_time} seconds")

    print("Conformer searching with force field done.")

    print("Making input files for semiempirical optimization...")

    semiempirical_opt_dir = os.path.join(output_dir, args.semiempirical_opt_folder)
    os.makedirs(semiempirical_opt_dir, exist_ok=True)
    inputs_dir = os.path.join(semiempirical_opt_dir, "inputs")
    outputs_dir = os.path.join(semiempirical_opt_dir, "outputs")
    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)

    for mol_id, smi in mol_ids_smis[args.task_id:len(mol_ids_smis):args.num_tasks]:
        ids = mol_id // 1000
        subinputs_dir = os.path.join(semiempirical_opt_dir, "inputs", f"inputs_{ids}")
        suboutputs_dir = os.path.join(semiempirical_opt_dir, "outputs", f"outputs_{ids}")
        os.makedirs(suboutputs_dir, exist_ok=True)
        if not os.path.exists(os.path.join(suboutputs_dir, f"{mol_id}.tar")) and os.path.exists(os.path.join(FF_conf_dir, "outputs", f"outputs_{ids}", f"{mol_id}_confs.sdf")):
            os.makedirs(subinputs_dir, exist_ok=True)
            if not os.path.exists(os.path.join(subinputs_dir, f"{mol_id}.in")) and not os.path.exists(os.path.join(subinputs_dir, f"{mol_id}.tmp")):
                with open(os.path.join(subinputs_dir, f"{mol_id}.in"), "w") as f:
                    f.write("")
                print(f"Making input file for semiempirical optimization for {mol_id}...")

    print("Optimizing conformers with semiempirical method...")

    for _ in range(1):
        for subinputs_folder in os.listdir(os.path.join(semiempirical_opt_dir, "inputs")):
            ids = subinputs_folder.split("_")[1]
            subinputs_dir = os.path.join(semiempirical_opt_dir, "inputs", f"inputs_{ids}")
            suboutputs_dir = os.path.join(semiempirical_opt_dir, "outputs", f"outputs_{ids}")
            for input_file in os.listdir(subinputs_dir):
                if ".in" in input_file:
                    mol_id = int(input_file.split(".in")[0])
                    try:
                        os.rename(os.path.join(subinputs_dir, input_file), os.path.join(subinputs_dir, f"{mol_id}.tmp"))
                    except:
                        continue
                    else:
                        ids = mol_id // 1000
                        smi = mol_id_to_smi[mol_id]
                        charge = mol_id_to_charge[mol_id]
                        mult = mol_id_to_mult[mol_id]
                        print(f"Optimizing conformers with semiempirical method for {mol_id} {smi}...")

                        tmp_mol_dir = os.path.join(suboutputs_dir, f"{mol_id}")
                        os.makedirs(tmp_mol_dir, exist_ok=True)

                        mol_confs_sdf = os.path.join(FF_conf_dir, "outputs", f"outputs_{ids}", f"{mol_id}_confs.sdf")
                        mols = RDKitMol.FromFile(mol_confs_sdf)
                        mol_id_to_FF_opted_xyz_dict = {}
                        mol_id_to_FF_opted_xyz_dict[mol_id] = {}
                        for conf_id, mol in enumerate(mols):
                            mol_id_to_FF_opted_xyz_dict[mol_id][conf_id] = mol.ToXYZ()
                        
                        start_time = time.time()
                        semiempirical_opt(mol_id, charge, mult, mol_id_to_FF_opted_xyz_dict, XTB_PATH, RDMC_PATH, G16_PATH, args.gaussian_semiempirical_opt_theory, args.gaussian_semiempirical_opt_n_procs, args.gaussian_semiempirical_opt_job_ram, args.scratch_dir, tmp_mol_dir, suboutputs_dir, subinputs_dir)
                        end_time = time.time()
                        print(f"Time for semiempirical optimization for {mol_id} took {end_time - start_time} seconds")

    print("Semiempirical optimization done.")

    print("Making input files for DFT optimization and frequency calculation")

    DFT_opt_freq_dir = os.path.join(output_dir, args.DFT_opt_freq_folder)
    os.makedirs(DFT_opt_freq_dir, exist_ok=True)
    inputs_dir = os.path.join(DFT_opt_freq_dir, "inputs")
    outputs_dir = os.path.join(DFT_opt_freq_dir, "outputs")
    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)

    for mol_id, smi in mol_ids_smis[args.task_id:len(mol_ids_smis):args.num_tasks]:
        ids = mol_id // 1000
        subinputs_dir = os.path.join(DFT_opt_freq_dir, "inputs", f"inputs_{ids}")
        suboutputs_dir = os.path.join(DFT_opt_freq_dir, "outputs", f"outputs_{ids}")
        os.makedirs(suboutputs_dir, exist_ok=True)
        semiempirical_opt_tar = os.path.join(semiempirical_opt_dir, "outputs", f"outputs_{ids}", f"{mol_id}.tar")
        if not os.path.exists(os.path.join(suboutputs_dir, f"{mol_id}", f"{mol_id}.log")) and os.path.exists(semiempirical_opt_tar):
            os.makedirs(subinputs_dir, exist_ok=True)
            if not os.path.exists(os.path.join(subinputs_dir, f"{mol_id}.in")) and not os.path.exists(os.path.join(subinputs_dir, f"{mol_id}.tmp")):
                with open(os.path.join(subinputs_dir, f"{mol_id}.in"), "w") as f:
                    f.write("")
                print(f"Making input file for DFT optimization and frequency calculation for {mol_id}...")

    print("Optimizing lowest energy semiempirical opted conformer with DFT method...")

    DFT_opt_freq_theories = [args.DFT_opt_freq_theory, args.DFT_opt_freq_theory_backup]

    for _ in range(1):
        for subinputs_folder in os.listdir(os.path.join(DFT_opt_freq_dir, "inputs")):
            ids = subinputs_folder.split("_")[1]
            subinputs_dir = os.path.join(DFT_opt_freq_dir, "inputs", f"inputs_{ids}")
            suboutputs_dir = os.path.join(DFT_opt_freq_dir, "outputs", f"outputs_{ids}")
            for input_file in os.listdir(subinputs_dir):
                if ".in" in input_file:
                    mol_id = int(input_file.split(".in")[0])
                    try:
                        os.rename(os.path.join(subinputs_dir, input_file), os.path.join(subinputs_dir, f"{mol_id}.tmp"))
                    except:
                        continue
                    else:
                        ids = mol_id // 1000
                        smi = mol_id_to_smi[mol_id]
                        charge = mol_id_to_charge[mol_id]
                        mult = mol_id_to_mult[mol_id]

                        output_mol_dir = os.path.join(DFT_opt_freq_dir, "outputs", f"outputs_{ids}", f"{mol_id}")
                        os.makedirs(output_mol_dir, exist_ok=True)

                        print(f"Optimizing lowest energy semiempirical opted conformer with DFT method for {mol_id} {smi}...")
                        
                        semiempirical_opt_tar = os.path.join(semiempirical_opt_dir, "outputs", f"outputs_{ids}", f"{mol_id}.tar")
                        failed_job, valid_job = semiempirical_opt_parser(mol_id, smi, semiempirical_opt_tar)

                        if valid_job:
                            mol_id_to_semiempirical_opted_xyz = get_mol_id_to_semiempirical_opted_xyz(valid_job)

                            converged = dft_scf_opt(mol_id, mol_id_to_semiempirical_opted_xyz[mol_id], G16_PATH, args.DFT_opt_freq_theory, args.DFT_opt_freq_n_procs, args.DFT_opt_freq_job_ram, charge, mult, args.scratch_dir, subinputs_dir, output_mol_dir)

                            if not converged:
                                print(f"DFT optimization for {mol_id} failed. Trying to optimize lowest energy FF opted conformer with DFT method...")
                                mol_confs_sdf = os.path.join(FF_conf_dir, "outputs", f"outputs_{ids}", f"{mol_id}_confs.sdf")
                                mols = RDKitMol.FromFile(mol_confs_sdf)
                                mol_id_to_FF_opted_xyz_dict = {}
                                mol_id_to_FF_opted_xyz_dict[mol_id] = {}
                                for conf_id, mol in enumerate(mols):
                                    if conf_id == 0:
                                        mol_id_to_FF_opted_xyz_dict[mol_id] = mol.ToXYZ()
                                        break
                            
                                converged = dft_scf_opt(mol_id, mol_id_to_FF_opted_xyz_dict[mol_id], G16_PATH, args.DFT_opt_freq_theory_backup, args.DFT_opt_freq_n_procs, args.DFT_opt_freq_job_ram, charge, mult, args.scratch_dir, subinputs_dir, output_mol_dir)

                        else:
                            print(f"All semiempirical opted conformers failed for {mol_id}")
                            print(failed_job)
                            
                            print("Trying to optimize lowest energy FF opted conformer with DFT method...")
                            mol_confs_sdf = os.path.join(FF_conf_dir, "outputs", f"outputs_{ids}", f"{mol_id}_confs.sdf")
                            mols = RDKitMol.FromFile(mol_confs_sdf)
                            mol_id_to_FF_opted_xyz_dict = {}
                            mol_id_to_FF_opted_xyz_dict[mol_id] = {}
                            for conf_id, mol in enumerate(mols):
                                if conf_id == 0:
                                    mol_id_to_FF_opted_xyz_dict[mol_id] = mol.ToXYZ()
                                    break
                            
                            converged = dft_scf_opt(mol_id, mol_id_to_FF_opted_xyz_dict[mol_id], G16_PATH, args.DFT_opt_freq_theory_backup, args.DFT_opt_freq_n_procs, args.DFT_opt_freq_job_ram, charge, mult, args.scratch_dir, subinputs_dir, output_mol_dir)

    print("DFT optimization and frequency calculation done.")

print("Done!")