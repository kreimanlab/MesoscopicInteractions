#!/bin/bash
arr=($(ls -S /Volumes/cuenap_ssd/h5 | grep .h5))
c=1
N_PERM=10000
METRIC='sg'
for i in "${arr[@]}"
do
    y=${i%.h5}
    s=${y##*/}
    d=`expr $c % 2 + 1`
    printf '> Start graph.py: %s, device: %s, metric: %s\n' $s $d $METRIC
    python3 graph.py $s $METRIC $d &
    printf '> Start perm.py: %s, device: %s, metric: %s, nperm: %s\n' $s $d $METRIC $N_PERM
    python3 perm.py $s $METRIC $N_PERM $d &
    c=`expr $c + 1` 
done
