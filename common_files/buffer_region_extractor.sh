#!/bin/bash

export MM_LATTICE="aftermm.lat"
export BUFFER_LATTICE="add.lat"
export TMP_G7="group_7_tmp.lat"
export TMP_G8="group_8_tmp.lat"
export TMP_G9="group_9_tmp.lat"
grep " 7 " $MM_LATTICE > $TMP_G7
grep " 8 " $MM_LATTICE > $TMP_G8
grep " 9 " $MM_LATTICE > $TMP_G9
NUM_BUFFER=`cat $TMP_G7 $TMP_G8 $TMP_G9 | wc -l`
echo $NUM_BUFFER

sed "/Atoms # atomic/q" $MM_LATTICE > $BUFFER_LATTICE

echo "" >> $BUFFER_LATTICE
cat $BUFFER_LATTICE

cat $TMP_G7 $TMP_G8 $TMP_G9 >> $BUFFER_LATTICE

sed '3d' $BUFFER_LATTICE > fuga.txt
sed "3i $NUM_BUFFER atoms" fuga.txt  > $BUFFER_LATTICE
rm fuga.txt
