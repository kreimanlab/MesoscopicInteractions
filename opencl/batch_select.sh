#!/bin/bash

# SUBJECTS
#    "m00001"  "m00018"  "m00025"  "m00033"  "m00044"  "m00053"  "m00061"  "m00083"
#    "m00003"  "m00019"  "m00026"  "m00035"  "m00045"  "m00055"  "m00068"  "m00084"
#    "m00004"  "m00021"  "m00027"  "m00037"  "m00047"  "m00056"  "m00071"  "m00095"
#    "m00005"  "m00022"  "m00028"  "m00038"  "m00048"  "m00058"  "m00073"  "m00096"
#    "m00006"  "m00023"  "m00030"  "m00039"  "m00049"  "m00059"  "m00075"  "m00097"
#    "m00017"  "m00024"  "m00032"  "m00043"  "m00052"  "m00060"  "m00079"  "m00100"
#    "m00107"  "m00122"  "m00124"  "mSu"

declare -a arr=(
#    "m00001"            "m00025"  "m00033"  "m00044"  "m00053"  "m00061"  "m00083"
#    "m00003"  "m00019"  "m00026"  "m00035"  "m00045"  "m00055"  "m00068"  "m00084"
#              "m00021"  "m00027"  "m00037"  "m00047"  "m00056"  "m00071"  "m00095"
#    "m00005"  "m00022"  "m00028"  "m00038"  "m00048"  "m00058"  "m00073"  "m00096"
#    "m00006"  "m00023"  "m00030"  "m00039"  "m00049"  "m00059"  "m00075"  "m00097"
#              "m00024"  "m00032"  "m00043"  "m00052"  "m00060"  "m00079"  "m00100"
#    "m00107"  "m00122"  "m00124"  "mSu"
    
#                                               "mSu"  "m00124"  "m00122"  "m00107"
#    "m00100"  "m00079"  "m00060"  "m00052"  "m00043"  "m00032"  "m00024"
#    "m00097"  "m00075"  "m00059"  "m00049"  "m00039"  "m00030"  "m00023"  "m00006"
#    "m00096"  "m00073"  "m00058"  "m00048"  "m00038"  "m00028"  "m00022"  "m00005"
#    "m00095"  "m00071"  "m00056"  "m00047"  "m00037"  "m00027"  "m00021"
#    "m00084"  "m00068"  "m00055"  "m00045"  "m00035"  "m00026"  "m00019"  "m00003"
#    "m00083"  "m00061"  "m00053"  "m00044"  "m00033"  "m00025"            "m00001"

    "m00005"
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


# --- Mar 22, 2018 ---

#    "sc"
#    "pc"

#    "m00001"  "m00018"  "m00025"  "m00033"  "m00044"  "m00053"  "m00061"  "m00083"
#    "m00003"  "m00019"  "m00026"  "m00035"  "m00045"  "m00055"  "m00068"  "m00084"
#    "m00004"  "m00021"  "m00027"  "m00037"  "m00047"  "m00056"  "m00071"  "m00095"
#    "m00005"  "m00022"  "m00028"  "m00038"  "m00048"  "m00058"  "m00073"  "m00096"
#    "m00006"  "m00023"  "m00030"  "m00039"  "m00049"  "m00059"  "m00075"  "m00097"
#    "m00017"  "m00024"  "m00032"  "m00043"  "m00052"  "m00060"  "m00079"  "m00100"
#    "m00107"  "m00122"  "m00124"


# --- Mar 22, 2018 ---

#    "s"
#    "p"
#    "sg"

#              "m00018"  "m00025"  "m00033"                                       
#    "m00003"                                          "m00055"                   
#    "m00004"  "m00021"  "m00027"                                                 
#    "m00005"  "m00022"  "m00028"                                                 
#                                            "m00049"                             
#    "m00017"                                                                     
#              "m00122"         


