#!/bin/bash

export MM_LATTICE="afterboxrelax.cuppt.lat"
export CU_PPT="cuppt_alone.lat"
export TMP_G10="group_10_tmp.lat"

echo "type 10 being extracted..."
grep " 10 " $MM_LATTICE > $TMP_G10
NUM_BUFFER=`cat $TMP_G10 | wc -l`
echo $NUM_BUFFER

sed "/Atoms # atomic/q" $MM_LATTICE > $CU_PPT

echo "" >> $CU_PPT
cat $CU_PPT

cat $TMP_G10 >> $CU_PPT

sed '3d' $CU_PPT > fuga.txt
sed "3i $NUM_BUFFER atoms" fuga.txt  > $CU_PPT
rm fuga.txt
