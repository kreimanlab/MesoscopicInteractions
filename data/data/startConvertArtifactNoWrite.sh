#!/bin/bash

echo "--- startConvert ---"
if [ -z "$1" ]
    then
        echo "(!) Usage: ./startConvert.sh m*****"
        exit
fi

# Copy files from NAS to local
#echo "[!] Working with directory:"
#du -hs /mnt/cuenap/data/txt/$1
#rsync -avP /mnt/cuenap/data/txt/$1 ./txt/
#du -hs /mnt/synology_2/data/$1
#rsync -avP /mnt/synology_2/data/$1 ./txt/

# Convert to H5eeg
echo "cd txt2h5eeg"
cd txt2h5eeg

echo "python3 main.py $1"
python3 main.py $1

# Remove original txt files
#rm -rf ../txt/$1/*

# Remove artifacts
cd ../h5eeg
./rmArtifactsNoWrite.sh $1

echo "<!> All done."
