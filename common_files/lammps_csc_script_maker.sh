#!/bin/bash
module load cray-python/3.8.5.1
source ~/my_env_dir_python3/bin/activate
python $MYWORK/common_files/yaml_to_lammps.py > lammps_script_modifier_temp.sh
source lammps_script_modifier_temp.sh
rm lammps_script_modifier_temp.sh
