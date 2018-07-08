#!/bin/bash

rm -f /tmp/job*.out

spath='/Users/steve/sim/plasticity/02.sensitivity/scripts'

for i in `seq 1 $1`;
do
	$spath/jobBackground.bash $2 $3 $4 $5
done 
