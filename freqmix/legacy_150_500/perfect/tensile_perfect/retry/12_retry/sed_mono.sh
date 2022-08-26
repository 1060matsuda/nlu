#!/bin/bash
cp sscurve.txt sscurve.csv
sed -i  -e 's/ /,/g' sscurve.csv
cp sstest.txt sstest.csv
sed -i  -e 's/ /,/g' sstest.csv
cp ssminimize.txt ssminimize.csv
sed -i  -e 's/ /,/g' ssminimize.csv
