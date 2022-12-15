#!/bin/bash
module load cray-python/3.8.5.1
source ~/my_env_dir_python3/bin/activate

for i in `seq 10`
do
	mkdir $i
done

for i in `seq 10`
do
	cp default/* $i/
done

for i in `seq 10`
do
	cd $i
	python $MYWORK/common_files/makebcclat_Eng.py > lat_maker.log &
	cd ..
done
