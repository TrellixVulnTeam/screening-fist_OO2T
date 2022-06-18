#!/bin/sh

# data processing script that runs data_proc2.py and outputs 
# ../../data/exp0$i-raw.csv files for each run 
# which contain traces and their metadata including compound mapping
# not normalized or anything

for i in $(seq 5); do
	workingdir="../../lab/0$i.0/ "
	echo $workingdir
	./data_proc2.py $workingdir > ../../data/exp0$i-raw.csv ; 
done
