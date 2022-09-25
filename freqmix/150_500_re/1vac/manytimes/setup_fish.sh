#!/bin/bash

for i in (seq 100)
	mkdir $i
end
for i in (seq 100)
	cp default/* $i/
end
for i in (seq 100)
	cd $i
	python $MYWORK/common_files/makebcclat_Eng.py > lat_maker.log
	cd ..
end
