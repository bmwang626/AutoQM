import os
import sys
import pickle as pkl
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
from radical_workflow.parser.cosmo_parser import cosmo_parser


def main(input_smiles_path, output_file_name, n_jobs, solvent_path, output_dir):
    submit_dir = os.getcwd()

    df_solute = pd.read_csv(input_smiles_path)
    mol_ids = list(df_solute.id)
    if "smiles" in df_solute.columns:
        mol_smis = list(df_solute.smiles)
    elif "smi" in df_solute.columns:
        mol_smis = list(df_solute.smi)
    elif "rxn_smi" in df_solute.columns:
        mol_smis = list(df_solute.rxn_smi)
    else:
        raise ValueError(f"No smiles column in input file. headers: {df_solute.columns}")

    mol_id_to_mol_smi = dict(zip(mol_ids, mol_smis))

    df_solvent = pd.read_csv(solvent_path)
    solvent_name_to_smi = dict(zip(df_solvent.cosmo_name, df_solvent.smiles))

    tar_file_paths = []
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if file.endswith(".tar"):
                tar_file_paths.append(os.path.join(root, file))

    out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(
        delayed(cosmo_parser)(tar_file_path) for tar_file_path in tqdm(tar_file_paths)
    )
    failed_mol_ids = [mol_ids[i] for i in range(len(mol_ids)) if out[i] is None]
    out = [x for x in out if x is not None]
    for each_data_lists in tqdm(out):
        for each_data_list in each_data_lists:
            for each_data in each_data_list:
                each_data[1] = solvent_name_to_smi[each_data[0]]
                each_data[3] = mol_id_to_mol_smi[each_data[2]]

    headers = [
        "solvent_name",
        "solvent_smiles",
        "solute_name",
        "solute_smiles",
        "temp (K)",
        "H (bar)",
        "ln(gamma)",
        "Pvap (bar)",
        "Gsolv (kcal/mol)",
        "Hsolv (kcal/mol)",
    ]

    cosmo_data_dict = {header: [] for header in headers}
    for each_data_lists in tqdm(out):
        for each_data_list in each_data_lists:
            for each_data in each_data_list:
                for i, header in enumerate(headers):
                    cosmo_data_dict[header].append(each_data[i])

    df_cosmo = pd.DataFrame(cosmo_data_dict)

    with open(os.path.join(submit_dir, f"{output_file_name}.pkl"), "wb") as outfile:
        pkl.dump(df_cosmo, outfile, protocol=pkl.HIGHEST_PROTOCOL)

    df_failed = df_solute[df_solute.id.isin(failed_mol_ids)]
    df_failed.to_csv(
        os.path.join(submit_dir, f"{output_file_name}_failed.csv"), index=False
    )

    print(f"Total number of jobs: {len(mol_ids)}")
    print(f"Failed mol ids: {len(failed_mol_ids)}")
    print(failed_mol_ids)


if __name__ == "__main__":
    input_smiles_path = sys.argv[1]
    output_file_name = sys.argv[2]
    n_jobs = int(sys.argv[3])
    solvent_path = sys.argv[4]
    cosmo_output_dir = sys.argv[5]

    main(input_smiles_path, output_file_name, n_jobs, solvent_path, cosmo_output_dir)
