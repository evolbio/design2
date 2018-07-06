#!/bin/sh

while [ 1 ]
do
	rsync fisher:/Users/steve/sim/plasticity/02.sensitivity/log/\*.log /Users/steve/sim/plasticity/02.sensitivity/log/.
	sleep 10
done