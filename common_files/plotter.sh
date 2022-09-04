#!/bin/bash

find ./*vac/ -maxdepth 0 -type d | tr "\n" "," > target_folders_plotter.csv