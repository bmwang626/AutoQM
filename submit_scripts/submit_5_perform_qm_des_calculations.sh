#!/bin/bash -l

# LLsub ./submit_5_perform_qm_des_calculations.sh [8,48,1] -q spot-xeon-p8
# squeue -p spot-xeon-p8
# watch LLloadSpot

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURM_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

conda activate rdmc_env
which python

scratch_dir=$TMPDIR/$USER/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-$LLSUB_RANK-$LLSUB_SIZE
mkdir -p $scratch_dir
echo $scratch_dir

# G16
export g16root=/home/gridsan/groups/RMG/Software/gaussian
export PATH=$g16root/g16/:$g16root/gv:$PATH
. $g16root/g16/bsd/g16.profile
export GAUSS_SCRDIR=$scratch_dir
chmod 750 $GAUSS_SCRDIR

# NBO
export nboroot=/home/gridsan/groups/RMG/Software/NBO7/
export nboram="10gb"
export PATH=$nboroot/nbo7/:$nboroot:$nboroot/nbo7/bin:$PATH

# AUTOQM
AUTOQM_PATH=/home/gridsan/groups/RMG/Software/autoqm
export PYTHONPATH=$AUTOQM_PATH:$PYTHONPATH

input_file=/home/gridsan/hwpang/qmdata_shared/qm_des_hwpang_shihcheng_oscar/input/quantum_green_species_data_24march12b_input.csv

python $AUTOQM_PATH/scripts/calculation/perform_qm_des_calculations.py \
    --input_file $input_file \
    --scratch_dir $scratch_dir \
    --task_id $LLSUB_RANK \
    --num_tasks $LLSUB_SIZE \
    --RDMC_path $RDMC_PATH \
    --G16_path $g16root/g16 \ 
    --title_card 

rm -rf $scratch_dir


