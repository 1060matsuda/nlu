#!/bin/bash
cp sscurve.txt sscurve.csv
sed -i  -e 's/ /,/g' sscurve.csv
cp outp2.txt outp2.csv
sed -i  -e 's/ /,/g' outp2.csv
