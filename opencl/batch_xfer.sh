#!/bin/bash

while :
do
    rsync -avP --remove-source-files ./results/*.h5 /Volumes/RawData/data/results/
    sleep 5
done
