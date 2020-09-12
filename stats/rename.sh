#!/bin/bash

declare -a arr=(
    "sub1"
    "sub2"
    "sub3"
    "sub4"
    "sub5"
    "sub6"
    "sub7"
    "sub8"
    "sub9"
    "sub10"
    "sub11"
    "sub12"
    "sub13"
    "sub14"
    "sub15"
    "sub16"
    "sub17"
    "sub18"
    "sub19"
    "sub20"
    "sub21"
    "sub22"
    "sub23"
    "sub24"
    "sub25"
    "sub26"
    "sub27"
    "sub28"
    "sub29"
    "sub30"
    "sub31"
    "sub32"
    "sub33"
    "sub34"
    "sub35"
    "sub36"
    "sub37"
    "sub38"
    "sub39"
    "sub40"
    "sub41"
    "sub42"
    "sub43"
    "sub44"
    "sub45"
    "sub46"
    "sub47"
    "sub48"
)

c=1
for sname in "${arr[@]}"
do
    sn="sub$c"
    echo $sn

    #cd $sn

    # rename all files
    #echo find . -iname "*$sname*" -exec rename "'s/$sname/$sn/'" '{}' \;
    #find . -iname "*$sname*" -exec rename "s/$sname/$sn/" '{}' \;
    #rm $sn

    # edit m0000* numbers to sub* in all files
    #echo find -name '*' -exec sed -i "'s/$sname/$sn/g'" {} \;
    find -name '*.m' -exec sed -i "s/$sname/$sn/g" '{}' \;
    find -name '*.txt' -exec sed -i "s/$sname/$sn/g" '{}' \;

    #cd ..

    c=`expr $c + 1`
done
