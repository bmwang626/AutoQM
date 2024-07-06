#!/bin/bash -l

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

conda activate autoqm_env
which python

#RDMC
RDMC_PATH=/home/gridsan/groups/RMG/Software/RDMC-main
export PATH=$RDMC_PATH:$PATH
export PYTHONPATH=$RDMC_PATH:$PYTHONPATH

#QMD
QMD_PATH=/home/gridsan/hwpang/qmdata_shared/quantum_green_reaction_complex_hwpang_shihcheng/AutoQM
export PYTHONPATH=$QMD_PATH:$PYTHONPATH

#input
input_smiles=$1

scratch_dir=$TMPDIR/$USER/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-$LLSUB_RANK-$LLSUB_SIZE
mkdir -p $scratch_dir
echo $scratch_dir

#r p complex semi opt
python $QMD_PATH/scripts/r_p_complex/reset_r_p_complex.py --input_smiles $input_smiles --RDMC_path $RDMC_PATH --scratch_dir $scratch_dir --task_id $LLSUB_RANK --num_tasks $LLSUB_SIZE

rm -rf $scratch_dir
