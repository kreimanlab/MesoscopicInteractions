#!/bin/bash

python3 clearNames.py $1
mkdir artifact_removed/$1
rsync -rvP --inplace $1/*.hdf5 artifact_removed/$1/
cd ../removeArtifacts/
python3 mainNoWrite.py $1
