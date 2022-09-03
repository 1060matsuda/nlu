#!/bin/bash
for i in `seq 100`
do
	cd $i
	bash sed_mono.sh
	cd ..
done

