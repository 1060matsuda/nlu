#!/bin/bash
cp outp_1.txt outp_1.csv
sed -i  -e 's/ /,/g' outp_1.csv
