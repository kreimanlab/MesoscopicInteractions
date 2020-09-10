#!/bin/bash


    #"*.sh"
declare -a arr=(
    "*.py"
    "*.m"
    "*.cl"
#    "*.sh"
)
for i in "${arr[@]}"
do
    #echo "> seldon - $i"
    #rsync -avP ./$i seldon:~/data/opencl/
    #echo "> leibniz - $i"
    #rsync -avP ./$i leibniz:~/data/opencl/
    echo "> cuebuntu - $i"
    rsync -avP ./$i cuebuntu:/nas_share/cuenap/scripts/opencl/
    echo "> cuebuntu2 - $i"
    rsync -avP ./$i cuebuntu2:/nas_share/RawData/scripts/opencl/
    #echo "> o2 - $i"
    #rsync -avP ./$i o2:/n/scratch2/jw324/opencl/
    #echo "> archimedes - $i"
    #rsync -avP ./$i archimedes:~/data/opencl/
done
