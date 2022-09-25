#!/bin/bash

nohup bash setup-run-analysis.sh &
echo $HOSTNAME > nohup_node.txt
