#!/bin/bash -l
#SBATCH -o slurm-orca_%A-%a.out
#SBATCH --ntasks=24
#SBATCH --array=0-31
#SBATCH --exclusive
#SBATCH --mem=0

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

echo $HOME

conda activate rdmc_env
which python

#RDMC
RDMC_PATH=/home/gridsan/groups/RMG/Software/RDMC-main
export PATH=$RDMC_PATH:$PATH
export PYTHONPATH=$PYTHONPATH:$RDMC_PATH

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/radical_workflow
export PYTHONPATH=$PYTHONPATH:$QMD_PATH

#openmpi
ompi=/home/gridsan/groups/RMG/Software/ompi-3.1.4
PATH=$ompi/bin:$PATH
LD_LIBRARY_PATH=$ompi/lib:$ompi/etc:$LD_LIBRARY_PATH

#Orca
orcadir=/home/gridsan/groups/RMG/Software/orca_4_2_1_linux_x86-64_openmpi314
export PATH=$PATH:$orcadir
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$orcadir
echo "orcaversion"
which orca
echo $PATH

#inputs
input_smiles=$1
echo $input_smiles
output_name=$2
echo $output_name
xyz_DFT_opt_dict=$3
echo $xyz_DFT_opt_dict
DLPNO_sp_n_procs=16
DLPNO_sp_job_ram=12000

echo "DLPNO_sp_n_procs $DLPNO_sp_n_procs"
echo "DLPNO_sp_job_ram $DLPNO_sp_job_ram"

#svp calculations
DLPNO_sp_folder="DLPNO_sp"
DLPNO_level_of_theory="uHF dlpno-ccsd(t) def2-svp def2-svp/c TightSCF NormalPNO"
echo $DLPNO_sp_folder
echo $DLPNO_level_of_theory

python -u $QMD_PATH/scripts/calculation/main_make_dlpno_sp_input.py \
--input_smiles $input_smiles \
--xyz_DFT_opt_dict $xyz_DFT_opt_dict \
--DLPNO_sp_folder $DLPNO_sp_folder \
--DLPNO_level_of_theory "$DLPNO_level_of_theory" \
--DLPNO_sp_n_procs $DLPNO_sp_n_procs \
--DLPNO_sp_job_ram $DLPNO_sp_job_ram \
--task_id $SLURM_ARRAY_TASK_ID \
--num_tasks $SLURM_ARRAY_TASK_COUNT \
--ORCA_path $orcadir

SubmitDir=`pwd`

echo "Starting DLPNO calculations..."

for i in {1..5}; do
    cd $SubmitDir
    for folder in $SubmitDir/output/$DLPNO_sp_folder/inputs/inputs_*; do

        folderind="${folder##$SubmitDir/output/$DLPNO_sp_folder/inputs/inputs_}"
        echo "folderind $folderind"
        
        for filename in $folder/*.in; do
            
            input=`basename $filename .in`
            
            if [ -e $folder/$input.in ]
            then
                echo "input $input"
                mv $folder/$input.in $folder/$input.tmp
                ScratchDir=$TMPDIR/$USER/orca/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-$input
                echo "ScratchDir $ScratchDir"
                mkdir -p $ScratchDir

                cd $ScratchDir
                echo $PATH
                cp $folder/$input.tmp $input.in
                $orcadir/orca $input.in > $input.log
                if [ -e $input.log ]
                then
                    if grep -Fq "ORCA TERMINATED NORMALLY" $input.log
                    then
                        echo "done"
                        cp $input.log $SubmitDir/output/$DLPNO_sp_folder/outputs/outputs_$folderind/
                        rm $folder/$input.tmp
                    elif grep -Fq "ORCA finished by error termination" $input.log
                    then
                        echo "done with error termination"
                        cp $input.log $SubmitDir/output/$DLPNO_sp_folder/outputs/outputs_$folderind/
                        rm $folder/$input.tmp
                    elif grep -Fq "The basis set was either not assigned or not available for this element" $input.log
                    then
                        echo "basis set not available"
                        cp $input.log $SubmitDir/output/$DLPNO_sp_folder/outputs/outputs_$folderind/
                        rm $folder/$input.tmp
                    elif grep -Fq "This wavefunction IS NOT FULLY CONVERGED!" $input.log
                    then
                        echo "wavefunction not converged"
                        cp $input.log $SubmitDir/output/$DLPNO_sp_folder/outputs/outputs_$folderind/
                        rm $folder/$input.tmp
                    else
                        echo "failed - unknown error"
                        cp $input.log $SubmitDir/output/$DLPNO_sp_folder/inputs/inputs_$folderind/
                        mv $folder/$input.tmp $folder/$input.in
                    fi
                else
                    echo "failed - no log file"
                    mv $folder/$input.tmp $folder/$input.in
                fi
                cd $SubmitDir
                rm -rf $ScratchDir
            fi
        done
    done
done

# f12 calculations
DLPNO_sp_folder="DLPNO_sp_f12"
DLPNO_level_of_theory="uHF UNO DLPNO-CCSD(T)-F12D cc-pvtz-f12 def2/J cc-pvqz/c cc-pvqz-f12-cabs RIJCOSX VeryTightSCF NormalPNO"
echo $DLPNO_sp_folder
echo $DLPNO_level_of_theory

python -u $QMD_PATH/scripts/calculation/main_make_dlpno_sp_input.py \
--input_smiles $input_smiles \
--xyz_DFT_opt_dict $xyz_DFT_opt_dict \
--DLPNO_sp_folder $DLPNO_sp_folder \
--DLPNO_level_of_theory "$DLPNO_level_of_theory" \
--DLPNO_sp_n_procs $DLPNO_sp_n_procs \
--DLPNO_sp_job_ram $DLPNO_sp_job_ram \
--ORCA_path $orcadir \
--task_id $SLURM_ARRAY_TASK_ID \
--num_tasks $SLURM_ARRAY_TASK_COUNT

SubmitDir=`pwd`

echo "Starting DLPNO calculations..."

for i in {1..5}; do
    cd $SubmitDir
    for folder in $SubmitDir/output/$DLPNO_sp_folder/inputs/inputs_*; do

        folderind="${folder##$SubmitDir/output/$DLPNO_sp_folder/inputs/inputs_}"
        echo "folderind $folderind"
        
        for filename in $folder/*.in; do
            
            input=`basename $filename .in`
            
            if [ -e $folder/$input.in ]
            then
                echo "input $input"
                mv $folder/$input.in $folder/$input.tmp
                ScratchDir=$TMPDIR/$USER/orca/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-$input
                echo "ScratchDir $ScratchDir"
                mkdir -p $ScratchDir

                cd $ScratchDir
                echo $PATH
                cp $folder/$input.tmp $input.in
                $orcadir/orca $input.in > $input.log
                if [ -e $input.log ]
                then
                    if grep -Fq "ORCA TERMINATED NORMALLY" $input.log
                    then
                        echo "done"
                        cp $input.log $SubmitDir/output/$DLPNO_sp_folder/outputs/outputs_$folderind/
                        rm $folder/$input.tmp
                    elif grep -Fq "ORCA finished by error termination" $input.log
                    then
                        echo "done with error termination"
                        cp $input.log $SubmitDir/output/$DLPNO_sp_folder/outputs/outputs_$folderind/
                        rm $folder/$input.tmp
                    elif grep -Fq "The basis set was either not assigned or not available for this element" $input.log
                    then
                        echo "basis set not available"
                        cp $input.log $SubmitDir/output/$DLPNO_sp_folder/outputs/outputs_$folderind/
                        rm $folder/$input.tmp
                    elif grep -Fq "This wavefunction IS NOT FULLY CONVERGED!" $input.log
                    then
                        echo "wavefunction not converged"
                        cp $input.log $SubmitDir/output/$DLPNO_sp_folder/outputs/outputs_$folderind/
                        rm $folder/$input.tmp
                    else
                        echo "failed - unknown error"
                        cp $input.log $SubmitDir/output/$DLPNO_sp_folder/inputs/inputs_$folderind/
                        mv $folder/$input.tmp $folder/$input.in
                    fi
                else
                    echo "failed - no log file"
                    mv $folder/$input.tmp $folder/$input.in
                fi
                cd $SubmitDir
                rm -rf $ScratchDir
            fi
        done
    done
done