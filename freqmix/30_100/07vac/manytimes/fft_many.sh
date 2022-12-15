#!/bin/bash
module load cray-python/3.8.5.1
source ~/my_env_dir_python3/bin/activate

echo "sed started"
for i in `seq 10`
do
	echo $i
	cd $i
	bash $MYWORK/common_files/sed_mono.sh > sed.log 
	cd ..
done
echo "FFT started"
for i in `seq 10`
do
	echo $i
	cd $i
	python $MYWORK/common_files/mono_detec_mix.py > FFT.log
	cd ..
done

echo "creating .gitignore"
cp $MYWORK/common_files/gitignore_betaonly ./default/.gitignore
for i in `seq 10`
do
	echo $i	
	cp ./default/.gitignore $i
done
