import sys
import pickle as pkl
import pandas as pd
from joblib import Parallel, delayed
import logging
import time
from tqdm import tqdm
import numpy as np

logging.basicConfig(level=logging.INFO)
start_time = time.time()

n_jobs = int(sys.argv[1])
logging.warning(f"n_jobs = {n_jobs}")
job_type = sys.argv[2]
only_cosmo = False
if len(sys.argv) == 4:
    if sys.argv[3] == "True":
        only_cosmo = True

if job_type == "reactants_products":
    projects = ["aug11b", "sep1a_filtered"]
elif job_type == "ts":
    projects = ["sep1a"]
logging.warning(f"projects = {projects}")

logging.warning("Loading inputs tables")
inputs_df_dict = {}
if job_type == "reactants_products":
    inputs_df_dict["aug11b"] = pd.read_csv(
        "./calculations/aug11b/inputs/reactants_products_aug11b_inputs.csv", index_col=0
    )
    inputs_df_dict["sep1a_filtered"] = pd.read_csv(
        "./calculations/sep1a_filtered/inputs/reactants_products_sep1a_filtered_inputs.csv",
        index_col=0,
    )
elif job_type == "ts":
    inputs_df_dict["sep1a"] = pd.read_csv(
        "./calculations/sep1a/inputs/ts_sep1a_inputs.csv", index_col=0
    )

mol_id_to_index_dict = {}
for project, df in inputs_df_dict.items():
    mol_id_to_index_dict[project] = {}
    for index, mol_id in zip(df.index, df.id):
        mol_id_to_index_dict[project][mol_id] = index

if not only_cosmo:

    for project, df in inputs_df_dict.items():
        df["project"] = [project for i in df.index]

    prop_names_to_remove = ["mol_smi"]

    def find_min_energy_conf(mol_dict):
        ens = np.array(
            [
                conf_dict["semiempirical_energy"]["scf"]
                for conf_id, conf_dict in mol_dict.items()
            ]
        )
        conf_ids = np.array([conf_id for conf_id, conf_dict in mol_dict.items()])
        lowest_conf_ind = conf_ids[np.argsort(ens)[0]]
        return lowest_conf_ind

    def fill_column(
        results_dict,
        inputs_df_dict,
        mol_id_to_index_dict,
        conf=False,
        semiempirical=False,
    ):
        for project, df in inputs_df_dict.items():

            success_dict = results_dict[project]

            prop_names = []
            columns_dict = {}
            for mol_id in success_dict:
                if conf:
                    for conf_id in success_dict[mol_id]:
                        for prop_name in success_dict[mol_id][conf_id]:
                            prop_names.append(prop_name)
                            for conf_id in range(10):
                                columns_dict[f"{prop_name}_conf_{conf_id}"] = [
                                    None for _ in df.index
                                ]
                else:
                    for prop_name in success_dict[mol_id]:
                        prop_names.append(prop_name)
                        columns_dict[prop_name] = [None for _ in df.index]
                    break
                break

            logging.warning("Creating columns")
            logging.warning(columns_dict.keys())

            if conf:
                for mol_id in tqdm(success_dict):
                    index = mol_id_to_index_dict[project][mol_id]
                    for conf_id in success_dict[mol_id]:
                        for prop_name in success_dict[mol_id][conf_id]:
                            columns_dict[f"{prop_name}_conf_{conf_id}"][index] = (
                                success_dict[mol_id][conf_id][prop_name]
                            )
            else:
                for mol_id in tqdm(success_dict):
                    index = mol_id_to_index_dict[project][mol_id]
                    for prop_name in success_dict[mol_id]:
                        columns_dict[prop_name][index] = success_dict[mol_id][prop_name]

            if semiempirical:
                out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(
                    delayed(find_min_energy_conf)(success_dict[mol_id])
                    for mol_id in success_dict
                )
                for prop_name in prop_names:
                    columns_dict[f"{prop_name}_min_energy_conf"] = [
                        None for _ in df.index
                    ]
                for mol_id, min_energy_conf_id in zip(success_dict.keys(), out):
                    index = mol_id_to_index_dict[project][mol_id]
                    for prop_name in success_dict[mol_id][min_energy_conf_id]:
                        columns_dict[f"{prop_name}_min_energy_conf"][index] = (
                            success_dict[mol_id][min_energy_conf_id][prop_name]
                        )

            df2 = pd.DataFrame(
                {
                    column_name: column
                    for column_name, column in columns_dict.items()
                    if column_name not in prop_names_to_remove
                }
            )
            df = pd.concat([df, df2], axis=1)
            inputs_df_dict[project] = df

        logging.warning("=" * 20)

    if job_type == "reactants_products":

        logging.warning("Loading ff results")
        ff_results_dict = {}
        with open(
            "./calculations/aug11b/reactants_products_aug11b_ff_opted_results.pkl", "rb"
        ) as f:
            ff_results_dict["aug11b"] = pkl.load(f)
        with open(
            "./calculations/sep1a_filtered/reactants_products_sep1a_filtered_ff_opted_results.pkl",
            "rb",
        ) as f:
            ff_results_dict["sep1a_filtered"] = pkl.load(f)

        logging.warning("Filling ff results")
        start_time_1 = time.time()
        fill_column(ff_results_dict, inputs_df_dict, mol_id_to_index_dict, conf=True)
        end_time_1 = time.time()
        logging.warning(f"Time taken: {end_time_1 - start_time_1}")

        logging.warning("Columns")
        logging.warning(inputs_df_dict[projects[0]].columns)

        logging.warning("Loading semiempirical results")
        semi_results_dict = {}
        with open(
            "./calculations/aug11b/reactants_products_aug11b_semiempirical_opted_results.pkl",
            "rb",
        ) as f:
            semi_results_dict["aug11b"] = pkl.load(f)
        with open(
            "./calculations/sep1a_filtered/reactants_products_sep1a_filtered_semiempirical_opted_results.pkl",
            "rb",
        ) as f:
            semi_results_dict["sep1a_filtered"] = pkl.load(f)

        logging.warning("Filling semiempirical results")
        start_time_1 = time.time()
        fill_column(
            semi_results_dict,
            inputs_df_dict,
            mol_id_to_index_dict,
            conf=True,
            semiempirical=True,
        )
        end_time_1 = time.time()
        logging.warning(f"Time taken: {end_time_1 - start_time_1}")

        logging.warning("Columns")
        logging.warning(inputs_df_dict[projects[0]].columns)

        logging.warning("Loading dft results")
        dft_results_dict = {}
        with open(
            "./calculations/aug11b/reactants_products_aug11b_dft_opted_results.pkl",
            "rb",
        ) as f:
            dft_results_dict["aug11b"] = pkl.load(f)
        with open(
            "./calculations/sep1a_filtered/reactants_products_sep1a_filtered_dft_opted_results.pkl",
            "rb",
        ) as f:
            dft_results_dict["sep1a_filtered"] = pkl.load(f)

        logging.warning("Filling dft results")
        start_time_1 = time.time()
        fill_column(dft_results_dict, inputs_df_dict, mol_id_to_index_dict, conf=False)
        end_time_1 = time.time()
        logging.warning(f"Time taken: {end_time_1 - start_time_1}")

        logging.warning("Columns")
        logging.warning(inputs_df_dict[projects[0]].columns)

    logging.warning("Loading dlpno results")
    dlpno_results_dict = {}
    if job_type == "reactants_products":
        with open(
            "./calculations/aug11b/reactants_products_aug11b_dlpno_sp_results.pkl", "rb"
        ) as f:
            dlpno_results_dict["aug11b"] = pkl.load(f)
        with open(
            "./calculations/sep1a_filtered/reactants_products_sep1a_filtered_dlpno_sp_results.pkl",
            "rb",
        ) as f:
            dlpno_results_dict["sep1a_filtered"] = pkl.load(f)
    elif job_type == "ts":
        with open("./calculations/sep1a/ts_sep1a_dlpno_sp_results.pkl", "rb") as f:
            dlpno_results_dict["sep1a"] = pkl.load(f)

    logging.warning("Filling dlpno results")
    start_time_1 = time.time()
    fill_column(dlpno_results_dict, inputs_df_dict, mol_id_to_index_dict, conf=False)
    end_time_1 = time.time()
    logging.warning(f"Time taken: {end_time_1 - start_time_1}")

    logging.warning("Columns")
    logging.warning(inputs_df_dict[projects[0]].columns)

    logging.warning("Concatenating hashed tables")
    results_df_merged = pd.concat(list(inputs_df_dict.values()), ignore_index=True)

    if job_type == "reactants_products":
        results_df_merged["dft_input_xyz_source"] = [
            None for _ in results_df_merged.index
        ]
        results_df_merged.loc[
            results_df_merged["dft_input_xyz"] != None, "dft_input_xyz_source"
        ] = "semiempirical_xyz_min_energy_conf"
        results_df_merged.loc[
            (results_df_merged["ff_xyz_conf_0"] == results_df_merged["dft_input_xyz"]),
            "dft_input_xyz_source",
        ] = "ff_xyz_conf_0"

    logging.warning("Saving all results table")
    if job_type == "reactants_products":
        with open(
            "./calculations/reactants_products_aug11b_sep1a_filtered_gfnff_xtb_wb97xd_dlpno_results_table.pkl",
            "wb",
        ) as f:
            pkl.dump(results_df_merged, f, protocol=pkl.HIGHEST_PROTOCOL)
    elif job_type == "ts":
        with open("./calculations/ts_sep1a_dlpno_results_table.pkl", "wb") as f:
            pkl.dump(results_df_merged, f, protocol=pkl.HIGHEST_PROTOCOL)

logging.warning("Loading cosmo results")
cosmo_results_dict = {}
if job_type == "reactants_products":
    with open(
        "./calculations/aug11b/reactants_products_aug11b_cosmo_results.pkl", "rb"
    ) as f:
        cosmo_results_dict["aug11b"] = pkl.load(f)
    with open(
        "./calculations/sep1a_filtered/reactants_products_sep1a_filtered_cosmo_results.pkl",
        "rb",
    ) as f:
        cosmo_results_dict["sep1a_filtered"] = pkl.load(f)
elif job_type == "ts":
    with open("./calculations/sep1a/ts_sep1a_cosmo_results.pkl", "rb") as f:
        cosmo_results_dict["sep1a"] = pkl.load(f)


def fill_column_cosmo(results_dict, hased_table_df_dict, mol_id_to_index_dict):
    props = [
        "H (bar)",
        "ln(gamma)",
        "Pvap (bar)",
        "Gsolv (kcal/mol)",
        "Hsolv (kcal/mol)",
    ]
    for project, df in inputs_df_dict.items():
        cosmo_result = results_dict[project]

        solvent_names = cosmo_result["solvent_name"].unique()
        temps = cosmo_result["temp (K)"].unique()

        columns_dict = {
            f"{solvent_name}_{temp}_{prop}": [None for _ in df.index]
            for solvent_name in solvent_names
            for temp in temps
            for prop in props
        }
        for row_ind in tqdm(cosmo_result.index):
            solvent_name = cosmo_result.loc[row_ind, "solvent_name"]
            temp = cosmo_result.loc[row_ind, "temp (K)"]
            mol_id = cosmo_result.loc[row_ind, "solute_name"]
            index = mol_id_to_index_dict[project][mol_id]
            for prop in props:
                columns_dict[f"{solvent_name}_{temp}_{prop}"][index] = cosmo_result.loc[
                    row_ind, prop
                ]

        for column_name, column in columns_dict.items():
            df[column_name] = column


def combine_cosmo_results(results_dict, mol_id_to_index_dict):
    header = [
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

    df = pd.DataFrame(columns=header)
    for project, cosmo_result in results_dict.items():
        cosmo_result["project"] = project
        df = pd.concat([df, cosmo_result], axis=0)
    return df


# logging.warning("Filling cosmo results")
# start_time_1 = time.time()
# fill_column_cosmo(cosmo_results_dict, inputs_df_dict, mol_id_to_index_dict)
# end_time_1 = time.time()
# logging.warning(f"Time taken: {end_time_1 - start_time_1}")
# logging.warning("Columns")
# logging.warning(inputs_df_dict[projects[0]].columns)

logging.warning("Combining cosmo results")
start_time_1 = time.time()
cosmo_df_merged = combine_cosmo_results(cosmo_results_dict, mol_id_to_index_dict)
end_time_1 = time.time()
logging.warning(f"Time taken: {end_time_1 - start_time_1}")

logging.warning("Saving cosmo results table")
if job_type == "reactants_products":
    with open(
        "./calculations/reactants_products_aug11b_sep1a_filtered_cosmo_results_table.pkl",
        "wb",
    ) as f:
        pkl.dump(cosmo_df_merged, f, protocol=pkl.HIGHEST_PROTOCOL)
elif job_type == "ts":
    with open("./calculations/ts_sep1a_cosmo_results_table.pkl", "wb") as f:
        pkl.dump(cosmo_df_merged, f, protocol=pkl.HIGHEST_PROTOCOL)

logging.warning("Done")
end_time = time.time()
logging.warning(f"Total time taken: {end_time - start_time}")
