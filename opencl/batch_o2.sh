#!/bin/bash

declare -a arr=(
    "sub1"  "sub10"  "sub16"  "sub22"  "sub28"  "sub34"  "sub40"
    "sub2"  "sub5"  "sub11"  "sub17"  "sub23"  "sub29"  "sub35"  "sub41"
    "sub6"  "sub12"  "sub18"  "sub24"  "sub30"  "sub36"  "sub42"
    "sub3"  "sub7"  "sub13"  "sub19"  "sub25"  "sub31"  "sub37"  "sub43"
    "sub4"  "sub8"  "sub14"  "sub20"  "sub26"  "sub32"  "sub38"  "sub44"
    "sub9"  "sub15"  "sub21"  "sub27"  "sub33"  "sub39"  "sub45"
    "sub46"  "sub47"  "sub48"  "mSu"
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
    "sd"
    "st"
    "sa"
    "sb"
    "sg"
)

# PERMUTATION PARAMS
N_PERM=10000
DO_PERM=true

# GRAPH PARAMS
DO_GRAPH=true

# SBATCH
TIME='119:00:00'
PYENV='/n/scratch2/jw324/opencl/pyenv1/bin/activate'
LOGDIR='/n/scratch2/jw324/opencl/log'
d='cpu'
MEM='8G'
p='medium'

for s in "${arr[@]}"
do
    c=$(( ( RANDOM % 2 ) ))
    for m in "${arr2[@]}"
    do
        #d=`expr $c % 2 + 1`
        printf '> Start graph.py: %s, device: %s, metric: %s\n' $s $d $m
        if [ $DO_GRAPH = true ] ; then
            printf "#!/bin/bash\n#SBATCH -o $LOGDIR/%%j-graph-$s-$m-$d.log\n#SBATCH -c 1\n#SBATCH -t $TIME\n#SBATCH -p $p\n#SBATCH --mem=$MEM\n\nsource $PYENV\npython3 graph.py $s $m $d\n" > cpujob.sh
            #sbatch cpujob.sh
        fi

        printf '> Start perm.py: %s, device: %s, metric: %s, nperm: %s\n' $s $d $m $N_PERM
        if [ $DO_PERM = true ] ; then
            printf "#!/bin/bash\n#SBATCH -o $LOGDIR/%%j-perm-$s-$m-$N_PERM-$d.log\n#SBATCH -c 1\n#SBATCH -t $TIME\n#SBATCH -p $p\n#SBATCH --mem=$MEM\n\nsource $PYENV\npython3 perm.py $s $m $N_PERM $d\n" > cpujob.sh
            #sbatch cpujob.sh
        fi
        c=`expr $c + 1` 
    done
done

echo "> All Done."


