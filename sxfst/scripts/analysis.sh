#!/bin/sh

./config.py -i ../../lab/01.0/ > ../../lab/01.0/config.yml
./config.py -i ../../lab/02.0/ > ../../lab/02.0/config.yml

for i in 3 4 5; do
	./config.py -i ../../lab/0$i.0/ -c ../../lab/01.0/config.yml  > ../../lab/0$i.0/config.yml
done

for i in $(seq 5); do
	workingdir="../../lab/0$i.0/config.yml "
	echo $workingdir
	./data_proc.py $workingdir
done
