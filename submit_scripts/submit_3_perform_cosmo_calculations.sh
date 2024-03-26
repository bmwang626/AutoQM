#!/bin/bash -l

# LLsub ./submit_3_perform_nonts_cosmo_calculations.sh [85,48,1] -q spot-xeon-p8
# squeue -p spot-xeon-p8
# watch LLloadSpot

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

conda activate autoqm_env
which python

#COSMO
TURBODIR=/home/gridsan/groups/RMG/Software/TmoleX19/TURBOMOLE
source $TURBODIR/Config_turbo_env
COSMOTHERM_PATH=/home/gridsan/groups/RMG/Software/COSMOtherm2023
COSMO_DATABASE_PATH=/home/gridsan/groups/RMG/COSMO_database/COSMObase2023

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/autoqm
export PYTHONPATH=$QMD_PATH:$PYTHONPATH

input_smiles=inputs/ts_nfho_dft_opt_freq_round_1_smi_for_rp_calc.csv
xyz_DFT_opt_dict=ts_nho_dft_opt_freq_round_1_rp_dft_opt_freq_results/ts_nho_dft_opt_freq_round_1_rp_xyz_input_ori.pkl

scratch_dir=$TMPDIR/$USER/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-$LLSUB_RANK-$LLSUB_SIZE
mkdir -p $scratch_dir
echo $scratch_dir

sleep $LLSUB_RANK

python -u $QMD_PATH/scripts/calculation/perform_cosmo_calculations.py \
    --input_smiles $input_smiles \
    --xyz_DFT_opt_dict $xyz_DFT_opt_dict \
    --scratch_dir $scratch_dir \
    --COSMO_input_pure_solvents $QMD_PATH/common_solvent_list_final.csv \
    --COSMOtherm_path $COSMOTHERM_PATH \
    --COSMO_database_path $COSMO_DATABASE_PATH \
    --task_id $LLSUB_RANK \
    --num_tasks $LLSUB_SIZE

rm -rf $scratch_dir


