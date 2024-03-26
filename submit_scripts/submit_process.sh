#!/bin/bash -l
#SBATCH -J process
#SBATCH -n 48
#SBATCH -o slurm-process.out

conda activate autoqm_env

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/autoqm
input_path=<path to the csv file>
output_name=<output name>

#FF conf xyzs
python -u $QMD_PATH/scripts/parsing/process_FF_conf_result_parallel.py $input_path ${output_name}_ff_opted_results 48

#semiempirical
python -u $QMD_PATH/scripts/parsing/process_semi_opt_result_parallel.py $input_path ${output_name}_semiempirical_opted_results 48

#dft
python -u $QMD_PATH/scripts/parsing/process_dft_opt_freq_result_parallel.py $input_path ${output_name}_dft_opted_results 48

#dlpno
python -u $QMD_PATH/scripts/parsing/process_dlpno_sp_result_parallel.py $input_path ${output_name}_dlpno_sp_results 48

#cosmo
python -u $QMD_PATH/scripts/parsing/process_cosmo_result_parallel.py $input_path ${output_name}_cosmo_results 48 $QMD_PATH/common_solvent_list_final.csv

# r p complexes ts
python -u $QMD_PATH/scripts/parsing/process_r_p_complex_ff_opt_result_parallel.py $input_path ${output_name}_r_p_complexes_ts_results 48
