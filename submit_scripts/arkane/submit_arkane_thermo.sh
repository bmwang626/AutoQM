#!/bin/bash -l

conda activate rmg_rdmc_env_20230623_v2

python /home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/arkane_processing/autoqm/scripts/arkane/calc_arkane_thermo.py \
    --csv_path /home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/arkane_processing/data/quantum_green_species_data_24march12b_for_arkane.csv \
    --energy_level QG-wb97xd/def2svp \
    --freq_level QG-wb97xd/def2svp \
    --energy_software gaussian \
    --freq_software gaussian \
    --freq_scale 0.9914
