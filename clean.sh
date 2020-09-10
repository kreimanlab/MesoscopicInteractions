#!/bin/bash


# Clear external files
read -p "Press [Enter] to clear files, including downloaded files..."
#echo "rm -rf ./data/h5eeg/example"
rm -rf ./data/h5eeg/*
rm -rf ./data/h5
rm -rf ./data/h5_notch20/example.h5
rm -rf ./data/h5_notch20/art_nosz/example_art.h5
rm -rf ./data/v2/verify
rm -rf ./data/v2/art
rm -rf ./data/v2/art_stim
rm -rf ./data/v2/art_nosz
rm -rf ./opencl/results/*
rm -rf ./stats/cache_debug/*
rm -rf ./stats/figures/*
echo "Done"


