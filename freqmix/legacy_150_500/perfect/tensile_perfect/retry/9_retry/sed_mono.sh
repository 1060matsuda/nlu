#!/bin/bash
cp sstest.txt sstest.csv
sed -i  -e 's/ /,/g' sstest.csv
cp ssminimize.txt ssminimize.csv
sed -i  -e 's/ /,/g' ssminimize.csv
