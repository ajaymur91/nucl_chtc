#!/bin/bash
echo 0 > /dev/shm/foo
for i in `seq 1 10`
do
	for j in `seq 1 10`
	do
		{ echo "$i - $j"; sleep 1; 
		if (( i==5 ))
		then
		echo 1 > /dev/shm/foo
		fi
	} &
	done
	wait
	if (( $(</dev/shm/foo) ))
	then
		echo "$i"
		break 2
	fi
done

echo hi
