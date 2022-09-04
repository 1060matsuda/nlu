#!/bin/bash

#Execute at the directory whose name is inf f1_f2 format, e.g., 20_100 or 150_500.
find ./*vac/ -maxdepth 0 -type d | tr "\n" "," > target_folders_plotter.csv
sed -i -e "s/,\$//" target_folders_plotter.csv
python $MYWORK/common_files/plotter_beta_vs_vac_rate.py
rm target_folders_plotter.csv