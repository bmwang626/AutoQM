#!/bin/bash -l

# LLsub ./submit_1_perform_nonts_ff_semi_dft_calculations.sh [1,3,16] -q spot-xeon-p8
# squeue -p spot-xeon-p8
# watch LLloadSpot

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

conda activate rdmc_env

#RDMC
RDMC_PATH=/home/gridsan/groups/RMG/Software/RDMC-main
export PATH=$RDMC_PATH:$PATH
export PYTHONPATH=$PYTHONPATH:$RDMC_PATH

#xtb
source /home/gridsan/groups/RMG/Software/xtb-6.4.1/share/xtb/config_env.bash
XTB_PATH=/home/gridsan/groups/RMG/Software/xtb-6.4.1/bin
export PATH=$XTB_PATH:$PATH

#gaussian
export PATH=$PATH:/home/gridsan/groups/RMG/Software/gaussian/g16
export PATH=$PATH:/home/gridsan/groups/RMG/Software/gaussian/gv
export g16root=/home/gridsan/groups/RMG/Software/gaussian
source /home/gridsan/groups/RMG/Software/gaussian/g16/bsd/g16.profile
export GAUSS_SCRDIR=""

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/autoqm
export PATH=$PATH:$QMD_PATH
export PYTHONPATH=$PYTHONPATH:$QMD_PATH

input_smiles=inputs/ts_nfho_dft_opt_freq_round_1_smi_for_rp_calc.csv

scratch_dir=$TMPDIR/$USER/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-$LLSUB_RANK-$LLSUB_SIZE
mkdir -p $scratch_dir
echo "Scratch directory: $scratch_dir"

python -u $QMD_PATH/scripts/calculation/perform_nonts_ff_semi_dft_calculations.py \
    --input_smiles $input_smiles \
    --XTB_path $XTB_PATH \
    --RDMC_path $RDMC_PATH \
    --G16_path $g16root/g16 \
    --scratch_dir $scratch_dir \
    --task_id $LLSUB_RANK \
    --num_tasks $LLSUB_SIZE

rm -rf $scratch_dir



