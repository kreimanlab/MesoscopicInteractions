#!/bin/bash

#declare -a arr=("m00006" "m00019" "m00023" "m00024" "m00026" "m00030" "m00037" "m00068" "m00043" "m00060" "m00083")
#declare -a arr=("m00019")
declare -a arr=("m00001" "m00003" "m00004" "m00005" "m00006" "m00017" "m00018" "m00019" "m00021" "m00022" "m00023" "m00024" "m00025" "m00026" "m00027" "m00028" "m00030" "m00032" "m00033" "m00035" "m00037" "m00038" "m00039" "m00043" "m00044" "m00045" "m00047" "m00048" "m00049" "m00052" "m00053" "m00055" "m00056" "m00058" "m00059" "m00060" "m00061" "m00068" "m00071" "m00073" "m00075" "m00079" "m00083" "m00084" "m00089" "m00095" "m00096" "m00097" "m00100" "m00107" "m00122" "m00124")

for i in "${arr[@]}"
do
    mkdir $i
    mkdir $i/surf
    mkdir $i/label
    rsync -IavP /mnt/cuenap_ssd/coregistration/$i/elec_recon ./$i/
    rsync -IavP /mnt/cuenap_ssd/coregistration/$i/surf/rh.pial ./$i/surf/
    rsync -IavP /mnt/cuenap_ssd/coregistration/$i/surf/lh.pial ./$i/surf/
    rsync -IavP /mnt/cuenap_ssd/coregistration/$i/label/*all_* ./$i/label/
    rsync -IavP /mnt/cuenap_ssd/coregistration/$i/label/*MACAQUE* ./$i/label/
    rsync -IavP /mnt/cuenap_ssd/coregistration/$i/label/*HCP* ./$i/label/
    rsync -IavP /mnt/cuenap_ssd/coregistration/$i/label/*.mat ./$i/label/
done
