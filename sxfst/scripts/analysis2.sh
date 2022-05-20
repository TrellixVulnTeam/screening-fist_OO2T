#!/bin/sh


for i in $(seq 5); do
	workingdir="../../lab/0$i.0/ "
	echo $workingdir
	./data_proc2.py $workingdir > ../../data/exp0$i-raw.csv ; 
done
