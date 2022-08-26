#!/bin/bash
cp outp_1.txt outp_1.csv
sed -i  -e 's/ /,/g' outp_1.csv
cp outp2.txt outp2.csv
sed -i  -e 's/ /,/g' outp2.csv
