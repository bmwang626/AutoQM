#!/bin/bash -l
#SBATCH -n 48
#SBATCH -N 1
#SBATCH --array=2-3

conda activate rmg_rdmc_env_20230623_v2

scratch_dir=$TMPDIR/$USER/$SLURM_JOB_ID/$SLURM_ARRAY_TASK_ID
mkdir -p $scratch_dir

if [ $SLURM_ARRAY_TASK_ID -eq 0 ]; then

    python /home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/arkane_processing/AutoQM/scripts/arkane/calc_arkane_thermo.py \
        --csv_path /home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/arkane_processing/data/quantum_green_species_data_24march12b_for_arkane.csv \
        --freq_level "qgwb97xd/def2svp" \
        --freq_software gaussian \
        --freq_scale 0.9914 \
        --energy_level "qgdlpnoccsd(t)f12d/ccpvtzf12" \
        --energy_software orca \
        --n_jobs 40 \
        --save_path /home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/arkane_processing/data/quantum_green_species_data_24march12b_dft_opted_dlpno_sp_thermos.csv \
        --scratch_dir $scratch_dir

elif [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then

    python /home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/arkane_processing/AutoQM/scripts/arkane/calc_arkane_thermo.py \
        --csv_path /home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/arkane_processing/data/quantum_green_species_data_24march12b_for_arkane.csv \
        --freq_level "qgwb97xd/def2svp" \
        --freq_software gaussian \
        --freq_scale 0.9914 \
        --energy_level "qgwb97xd/def2svp" \
        --energy_software gaussian \
        --n_jobs 40 \
        --save_path /home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/arkane_processing/data/quantum_green_species_data_24march12b_dft_opted_dft_sp_thermos.csv \
        --scratch_dir $scratch_dir

elif [ $SLURM_ARRAY_TASK_ID -eq 2 ]; then

    python /home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/arkane_processing/AutoQM/scripts/arkane/calc_arkane_thermo.py \
        --csv_path /home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/arkane_processing/data/quantum_green_species_data_24march12b_for_arkane.csv \
        --freq_level "qgwb97xd/def2svp" \
        --freq_software gaussian \
        --freq_scale 0.9914 \
        --energy_level "qgdlpnoccsd(t)f12d/ccpvtzf12" \
        --energy_software orca \
        --n_jobs 40 \
        --save_path /home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/arkane_processing/data/quantum_green_species_data_24march12b_dft_opted_dlpno_sp_thermos_no_bac.csv \
        --scratch_dir $scratch_dir \
        --no_bac_for_thermo

elif [ $SLURM_ARRAY_TASK_ID -eq 3 ]; then

    python /home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/arkane_processing/AutoQM/scripts/arkane/calc_arkane_thermo.py \
        --csv_path /home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/arkane_processing/data/quantum_green_species_data_24march12b_for_arkane.csv \
        --freq_level "qgwb97xd/def2svp" \
        --freq_software gaussian \
        --freq_scale 0.9914 \
        --energy_level "qgwb97xd/def2svp" \
        --energy_software gaussian \
        --n_jobs 40 \
        --save_path /home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/arkane_processing/data/quantum_green_species_data_24march12b_dft_opted_dft_sp_thermos_no_bac.csv \
        --scratch_dir $scratch_dir \
        --no_bac_for_thermo

fi

