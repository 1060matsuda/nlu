#!/bin/bash

module load python

for i in `seq 100`
do
cd $i
python $COMMONFILES/mono_detec_mix.py
cd ..
done
