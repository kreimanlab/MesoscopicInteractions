#!/bin/bash

declare -a arr=(
#    "sub1"            "sub10"  "sub16"  "sub22"  "sub28"  "sub34"  "sub40"
#    "sub2"  "sub5"  "sub11"  "sub17"  "sub23"  "sub29"  "sub35"  "sub41"
#              "sub6"  "sub12"  "sub18"  "sub24"  "sub30"  "sub36"  "sub42"
#    "sub3"  "sub7"  "sub13"  "sub19"  "sub25"  "sub31"  "sub37"  "sub43"
#    "sub4"  "sub8"  "sub14"  "sub20"  "sub26"  "sub32"  "sub38"  "sub44"
#              "sub9"  "sub15"  "sub21"  "sub27"  "sub33"  "sub39"  "sub45"
#    "sub46"  "sub47"  "sub48"  "mSu"
    
#                                               "mSu"  "sub48"  "sub47"  "sub46"
#    "sub45"  "sub39"  "sub33"  "sub27"  "sub21"  "sub15"  "sub9"
#    "sub44"  "sub38"  "sub32"  "sub26"  "sub20"  "sub14"  "sub8"  "sub4"
#    "sub43"  "sub37"  "sub31"  "sub25"  "sub19"  "sub13"  "sub7"  "sub3"
#    "sub42"  "sub36"  "sub30"  "sub24"  "sub18"  "sub12"  "sub6"
#    "sub41"  "sub35"  "sub29"  "sub23"  "sub17"  "sub11"  "sub5"  "sub2"
#    "sub40"  "sub34"  "sub28"  "sub22"  "sub16"  "sub10"            "sub1"

    "sub3"
)

# METRICS
# s     Spearman correlation coefficient
# p     Pearson correlation coefficient
# sP    Spearman partial correlation coefficient
# pP    Pearson partial correlation coefficient
# sc    Spearman magnitude squared coherence
# pc    Pearson magnitude squared coherence
# sd    Spearman correlation coefficient of delta envelope
# st    Spearman correlation coefficient of theta envelope
# sa    Spearman correlation coefficient of alpha envelope
# sb    Spearman correlation coefficient of beta envelope
# sg    Spearman correlation coefficient of gamma envelope
# shg    Spearman correlation coefficient of high gamma envelope
declare -a arr2=(
    "pc"
)

# PERMUTATION PARAMS
N_PERM=10000
DO_PERM=true

# GRAPH PARAMS
DO_GRAPH=false


for s in "${arr[@]}"
do
    c=$(( ( RANDOM % 2 ) ))
    for m in "${arr2[@]}"
    do
        #d=`expr $c % 2 + 1`
        d=0
        printf '> Start graph.py: %s, device: %s, metric: %s\n' $s $d $m
        if [ $DO_GRAPH = true ] ; then
            python3 graph.py $s $m $d
            #python3 graph.py $s $m $d &
        fi

        printf '> Start perm.py: %s, device: %s, metric: %s, nperm: %s\n' $s $d $m $N_PERM
        if [ $DO_PERM = true ] ; then
            python3 perm.py $s $m $N_PERM $d
            #python3 perm.py $s $m $N_PERM $d &
        fi
        c=`expr $c + 1` 
    done
done

echo "> All Done."


