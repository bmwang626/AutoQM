import os
from argparse import ArgumentParser

import pandas as pd
from joblib import Parallel, delayed
from rdkit import Chem
from rdmc.mol import RDKitMol


def sdf_parser(sdf_file_path):
    rxn_id = int(os.path.basename(sdf_file_path).split(".")[0].split("_")[1])

    r_mol, p_mol, ts_mol = RDKitMol.FromFile(
        sdf_file_path, removeHs=False, sanitize=False
    )

    r_smi = r_mol.GetProp("_Name")
    p_smi = p_mol.GetProp("_Name")
    rxn_smi = ts_mol.GetProp("_Name")

    pre_r_mol = RDKitMol.FromSmiles(r_smi, removeHs=False, sanitize=True)
    pre_p_mol = RDKitMol.FromSmiles(p_smi, removeHs=False, sanitize=True)

    pre_r_adj = pre_r_mol.GetAdjacencyMatrix()
    pre_p_adj = pre_p_mol.GetAdjacencyMatrix()

    r_adj = r_mol.GetAdjacencyMatrix()
    if not (r_adj == pre_r_adj).all():
        return None
    p_adj = p_mol.GetAdjacencyMatrix()
    if not (p_adj == pre_p_adj).all():
        return None

    return rxn_id, r_smi, p_smi, rxn_smi, r_mol._mol, p_mol._mol, ts_mol._mol

parser = ArgumentParser()
parser.add_argument(
    "--input_smiles",
    type=str,
    required=True,
    help="input smiles included in a .csv file",
)
parser.add_argument(
    "--output_name", type=str, required=True, help="output file name"
)
parser.add_argument(
    "--num_tasks",
    type=int,
    required=True,
)

args = parser.parse_args()
input_smiles_path = args.input_smiles
output_file_name = args.output_name
n_jobs = args.num_tasks

df = pd.read_csv(input_smiles_path)

rxn_ids = df.id

success_sdf_paths = []
for root, dirs, files in os.walk("output/r_p_complex_ff_opt/outputs"):
    for file in files:
        if file.endswith(".sdf"):
            success_sdf_paths.append(os.path.join(root, file))

out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(
    delayed(sdf_parser)(sdf_path) for sdf_path in success_sdf_paths
)
out = [x for x in out if x is not None]

print(f"Total number of reactions: {len(rxn_ids)}")
print(f"Number of successful TSs: {len(success_sdf_paths)}")
print(f"Number of successful r complexes, p complexes, TSs: {len(out)}")

success_rxn_ids = [x[0] for x in out]
success_r_smis = [x[1] for x in out]
success_p_smis = [x[2] for x in out]

csv_file = os.path.join(f"{output_file_name}.csv")
df_out = pd.DataFrame(
    {"id": success_rxn_ids, "rsmi": success_r_smis, "psmi": success_p_smis}
)
df_out.to_csv(csv_file, index=False)

r_writer = Chem.rdmolfiles.SDWriter(os.path.join(f"{output_file_name}_reactants.sdf"))
p_writer = Chem.rdmolfiles.SDWriter(os.path.join(f"{output_file_name}_products.sdf"))
ts_writer = Chem.rdmolfiles.SDWriter(os.path.join(f"{output_file_name}_ts.sdf"))

for rxn_id, r_smi, p_smi, rxn_smi, r_mol, p_mol, ts_mol in out:
    r_mol.SetProp("_Name", f"{rxn_id}_{r_smi}")
    p_mol.SetProp("_Name", f"{rxn_id}_{p_smi}")
    ts_mol.SetProp("_Name", f"{rxn_id}_{rxn_smi}")
    r_writer.write(r_mol)
    p_writer.write(p_mol)
    ts_writer.write(ts_mol)

r_writer.close()
p_writer.close()
ts_writer.close()
