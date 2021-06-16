#!/bin/bash
rm -rf batch.dag
for i in `seq 0 99`
do
cat << EOF >> batch.dag
JOB $i job.sub
VARS $i index="$i"
EOF
done
