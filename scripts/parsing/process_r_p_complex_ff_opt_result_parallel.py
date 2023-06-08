import os
import sys
import csv
import tarfile
import pickle as pkl
import pandas as pd
from joblib import Parallel, delayed

from rdkit import Chem
from rdmc.mol import RDKitMol


def parser(sdf_file_path):
    rxn_id = int(os.path.basename(sdf_file_path).split(".")[0].split("_")[1])

    r_mol, p_mol, ts_mol = RDKitMol.FromFile(
        sdf_file_path, removeHs=False, sanitize=False
    )

    r_smi = r_mol.GetProp("_Name")
    p_smi = p_mol.GetProp("_Name")

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

    return rxn_id, r_smi, p_smi, ts_mol._mol, r_mol._mol, p_mol._mol


input_smiles_path = sys.argv[1]
output_file_name = sys.argv[2]
n_jobs = int(sys.argv[3])

df = pd.read_csv(input_smiles_path)

rxn_ids = df.id

if "rxn_smi" in df.columns:
    rxn_smis = df.rxn_smi
elif "rxn_smiles" in df.columns:
    rxn_smis = df.rxn_smiles
else:
    raise ValueError("No reaction smiles provided")

success_sdf_paths = []
for root, dirs, files in os.walk("output/r_p_complex_ff_opt/outputs"):
    for file in files:
        if file.endswith(".sdf"):
            success_sdf_paths.append(os.path.join(root, file))

out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(
    delayed(parser)(sdf_path) for sdf_path in success_sdf_paths
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

ts_writer = Chem.rdmolfiles.SDWriter(os.path.join(f"{output_file_name}_ts.sdf"))
r_writer = Chem.rdmolfiles.SDWriter(os.path.join(f"{output_file_name}_reactants.sdf"))
p_writer = Chem.rdmolfiles.SDWriter(os.path.join(f"{output_file_name}_products.sdf"))

for rxn_id, r_smi, p_smi, ts_mol, r_mol, p_mol in out:
    ts_writer.write(ts_mol)
    r_writer.write(r_mol)
    p_writer.write(p_mol)

ts_writer.close()
r_writer.close()
p_writer.close()
