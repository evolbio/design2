#!/bin/bash

num=$(( ( RANDOM % 10000 )  + 1 ))

nohup $1 $2 $3 $4 > /tmp/jobout$num.out 2> /tmp/joberror$num.out < /dev/null &
