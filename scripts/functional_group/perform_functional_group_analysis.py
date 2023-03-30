import os
from tqdm import tqdm
import pandas as pd
from argparse import ArgumentParser
from joblib import Parallel, delayed

from radical_workflow.functional_group.analysis import functional_group_analysis

parser = ArgumentParser()
parser.add_argument('--input_smiles', type=str, required=True,
                    help='input smiles included in a .csv file')
parser.add_argument('--n_jobs', type=int, default=1,
                    help='number of workers to use')

args = parser.parse_args()

# load smiles
print('Loading smiles...')
df_smiles = pd.read_csv(args.input_smiles, index_col=0)
smiles_list = df_smiles['smiles'].tolist()

# perform functional group analysis
print('Performing functional group analysis...')
outs = Parallel(n_jobs=args.n_jobs, backend="multiprocessing")(delayed(functional_group_analysis)(smi) for smi in tqdm(smiles_list))

# collect results
functional_group_smiles = set()
for out in outs:
    functional_group_smiles.update(out)

# save results
print('Saving results...')
df_groups = pd.DataFrame(functional_group_smiles, columns=['functional_group_smiles'])
input_file_name = os.path.basename(args.input_smiles)
output_file_name = input_file_name.replace('.csv', '_functional_groups.csv')
df_groups.to_csv(output_file_name)