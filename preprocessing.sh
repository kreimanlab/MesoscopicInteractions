#!/bin/bash
echo source downloads/pyenv/bin/activate
source downloads/pyenv/bin/activate

echo ""
echo "  Preprocessing"
echo ""
PS3=": "

select opt in "run txt2h5eeg" "run artifact removal" "run synth" "run artifact removal (v2)" quit; do

  case $opt in
    "run txt2h5eeg")

echo ""
echo "      This script converts exported .txt files from Natus Neuroworks 8"
echo "      software (neuro.natus.com). The short export file located at"
echo "      ./data/txt/example/_Export-example_001.txt will be converted to an"
echo "      HDF5-formatted file. A notch filter is applied to remove power line"
echo "      noise, using the script mainfilter.m in the folder ./data/txt2h5eeg"
echo ""

      echo "[*] This script should be configured to run with the following directories:"
      echo "In ./data/txt2h5eeg/config:"
      echo "    EXPORTED_DIR=../txt"
      echo "    FINISHED_DIR=../h5eeg"
      echo ""
      read -p "Press [Enter] to run."
      currDir=`pwd`
      echo "cd data/txt2h5eeg"
      cd "data/txt2h5eeg"
      echo "./main.py example"
      ./main.py example

      # return to starting dir
      echo "cd $currDir"
      cd $currDir

      echo "[*] The output is in: ./data/h5eeg/example/example_001.hdf5"
      break
      ;;

    "run artifact removal")
      echo "artifact removal"

echo ""
echo "      This script marks artifacts at 1-second intervals using the scripts"
echo "      in ./data/scripts/removeArtifacts"
echo ""

      currDir=`pwd`
      echo "mkdir ./data/h5eeg/artifact_removed"
      mkdir ./data/h5eeg/artifact_removed
      echo "cp data/txt2h5eeg/rmArtifactsNoWrite.sh data/h5eeg/"
      cp data/removeArtifacts/rmArtifactsNoWrite.sh data/h5eeg/
      echo "cp data/txt2h5eeg/clearNames.py data/h5eeg/"
      cp data/removeArtifacts/clearNames.py data/h5eeg/
      echo "cd ./data/h5eeg"
      cd ./data/h5eeg
      echo "./rmArtifactsNoWrite.sh example"
      ./rmArtifactsNoWrite.sh example

      # return to starting dir
      echo "cd $currDir"
      cd $currDir

      echo "[*] The output is in: ./data/h5eeg/artifact_removed/example"
      break
      ;;

    "run synth")
      echo "synth"
      currDir=`pwd`
      echo "mkdir ./data/h5"
      mkdir ./data/h5
      echo "cd ./data/synth"
      cd ./data/synth
      echo matlab -nodesktop -nosplash -r "synth;exit()"
      matlab -nodesktop -nosplash -r "synth;exit()"
      
      # return to starting dir
      echo "cd $currDir"
      cd $currDir

      break
      #read -p "Enter the first number: " n1
      ;;

    "run artifact removal (v2)")
      echo "artifact removal (v2)"
      currDir=`pwd`
      echo "cd ./data/v2"
      cd ./data/v2
      echo matlab -nodesktop -nosplash -r "notch20;annotate_fn('example');exit()"
      matlab -nodesktop -nosplash -r "notch20;annotate_fn('example');exit()"

      # make artifact annotation folder
      echo mkdir $currDir/data/h5_notch20/art
      mkdir $currDir/data/h5_notch20/art
      echo mkdir $currDir/data/v2/art
      mkdir $currDir/data/v2/art
      echo mkdir $currDir/data/v2/art_stim
      mkdir $currDir/data/v2/art_stim
      echo mkdir $currDir/data/v2/art_nosz
      mkdir $currDir/data/v2/art_nosz
      echo mv $currDir/data/h5_notch20/*_art.h5 $currDir/data/h5_notch20/art/
      mv $currDir/data/h5_notch20/*_art.h5 $currDir/data/h5_notch20/art/

      # matlab commands
      #echo matlab -nodesktop -nosplash -r "verify_final;exit()"
      #matlab -nodesktop -nosplash -r "verify_final;exit()"
      #echo matlab -nodesktop -nosplash -r "plot_sz;exit()"
      #matlab -nodesktop -nosplash -r "plot_sz;exit()"
      #echo matlab -nodesktop -nosplash -r "plot_sz2new;exit()"
      #matlab -nodesktop -nosplash -r "plot_sz2new;exit()"
      #echo matlab -nodesktop -nosplash -r "plot_sz2stim;exit()"
      #matlab -nodesktop -nosplash -r "plot_sz2stim;exit()"
      #echo matlab -nodesktop -nosplash -r "manual;exit()"
      #matlab -nodesktop -nosplash -r "manual;exit()"
      #echo matlab -nodesktop -nosplash -r "manual_stim;exit()"
      #matlab -nodesktop -nosplash -r "manual_stim;exit()"
      #echo matlab -nodesktop -nosplash -r "manual_sz;exit()"
      #matlab -nodesktop -nosplash -r "manual_sz;exit()"

      echo matlab -nodesktop -nosplash -r "verify_final;plot_sz;plot_sz2new;plot_sz2stim;manual;manual_stim;manual_sz;exit()"
      matlab -nodesktop -nosplash -r "verify_final;plot_sz;plot_sz2new;plot_sz2stim;manual;manual_stim;manual_sz;exit()"

      # return to starting dir
      echo "cd $currDir"
      cd $currDir

      # attach ground-truth neighboring electrode definitions
      echo cd ./data/mgrid
      cd ./data/mgrid
      echo matlab -nodesktop -nosplash -r "h5_add_mgrid;exit()"
      matlab -nodesktop -nosplash -r "h5_add_mgrid;exit()"
      echo "Finished h5_add_mgrid."
      
      # return to starting dir
      echo "cd $currDir"
      cd $currDir

      # copy artifact directory to finished .h5 directory
      echo rm -rf data/h5_notch20/art
      rm -rf data/h5_notch20/art
      echo rsync -avP data/v2/art_nosz data/h5_notch20/
      rsync -avP data/v2/art_nosz data/h5_notch20/
      
      # exit message
      echo "Done, artifact annotation summary figures are in: ./data/v2/verify/"
      
      break
      ;;

    quit)
      break
      ;;

    *) 
      echo "Invalid option $REPLY"
      ;;

  esac
done

