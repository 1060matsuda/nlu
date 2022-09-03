#!/bin/bash
read -p "This script will delete all the directories. Are you sure? (y/N): " yn
case "$yn" in 
[yY]*)
for i in `seq 100`
do
	rm -r $i
done
;; 
*) 
echo "aborting" 
;; 
esac
