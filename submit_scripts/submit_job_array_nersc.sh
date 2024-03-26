#!/bin/bash -l
#SBATCH -q regular
##SBATCH -q debug
#SBATCH -C haswell
#SBATCH -J cosmo
#SBATCH -L SCRATCH
#SBATCH -c 4
##SBATCH --array=0-1
#SBATCH -t 2:00:00
#SBATCH --array=0-199

source ~/.bashrc
conda activate QM_descriptors

#COSMO
#TURBODIR=/project/projectdirs/m1839/Software/COSMO/TmoleX2021/TURBOMOLE
TURBODIR=/project/projectdirs/m1839/Software/COSMO/TURBOMOLE
source $TURBODIR/Config_turbo_env
COSMOTHERM_PATH=/project/projectdirs/m1839/Software/COSMO/COSMOtherm2021
COSMO_DATABASE_PATH=/project/projectdirs/m1839/Software/COSMO/COSMObase2021

echo "My task ID: " $SLURM_ARRAY_TASK_ID
echo "Number of Tasks: " $SLURM_ARRAY_TASK_COUNT

#QMD
QMD_PATH=~/Software/QM_descriptors_calculation

input_smiles=descend_300K_data_down_30K_tree-20220520_sampled_mols.csv
xyz_DFT_opt=CNO_radicals_dft_opt_xyz.pkl

python $QMD_PATH/main.py --input_smiles $input_smiles --task_id $SLURM_ARRAY_TASK_ID --num_tasks $SLURM_ARRAY_TASK_COUNT --COSMOtherm_path $COSMOTHERM_PATH --COSMO_database_path $COSMO_DATABASE_PATH --skip_conf_search_FF --skip_semiempirical_opt --skip_DFT_opt_freq --xyz_DFT_opt $xyz_DFT_opt --COSMO_input_pure_solvents $QMD_PATH/common_solvent_list_final.csv --skip_DLPNO
# python $QMD_PATH/main.py --input_smiles $input_smiles --task_id $SLURM_ARRAY_TASK_ID --num_tasks $SLURM_ARRAY_TASK_COUNT --COSMOtherm_path $COSMOTHERM_PATH --COSMO_database_path $COSMO_DATABASE_PATH --skip_conf_search_FF --skip_semiempirical_opt --skip_DFT_opt_freq --skip_COSMO --xyz_DFT_opt $xyz_DFT_opt
~                                                                                                                                                                                                                                              
~                                                                                                                                                                                                                                              
~                                                                                                                                                                                                                                              
~                                                                                                                                                                                                                                              
~                                                                                                                                                                                                                                              
~                                                                                                                                                                                                                                              
~                                                                                                                                                                                                                                              
~                                                                                                                                                                                                                                              
~                                                                                                                                                                                                                                              
~          