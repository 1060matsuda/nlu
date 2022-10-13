#!/bin/bash
python $MYWORK/common_files/yaml_to_lammps_nrb_mix.py > lammps_script_modifier_temp.sh
source lammps_script_modifier_temp.sh
rm lammps_script_modifier_temp.sh
