#!/bin/bash
python $MYWORK/common_files/yaml_to_lammps.py > lammps_script_modifier_temp.sh
source lammps_script_modifier_temp.sh
rm lammps_script_modifier_temp.sh
