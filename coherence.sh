#!/bin/bash
echo source downloads/pyenv/bin/activate
source downloads/pyenv/bin/activate

echo ""
echo "  Coherence"
echo ""

PS3=": "

select opt in "run permutation" "run coherence" "run combine" quit; do

  case $opt in
    "run permutation")
      echo "permutation"
      currDir=`pwd`
      echo cd opencl
      cd opencl
      echo python3 perm.py example pc 1000 0
      python3 perm.py example pc 1000 0
      echo "Note: if failed, make sure pyopencl is properly installed and configured"
      # return to starting dir
      echo "cd $currDir"
      cd $currDir
      break
      ;;
    "run coherence")
      echo "coherence"
      currDir=`pwd`
      echo cd opencl
      cd opencl
      echo python3 graph.py example pc 0
      python3 graph.py example pc 0
      echo "Note: if failed, make sure pyopencl is properly installed and configured"
      # return to starting dir
      echo "cd $currDir"
      cd $currDir
      break
      ;;
    "run combine")
      echo "combine"
      currDir=`pwd`
      echo cd opencl
      cd opencl
      echo rm results/*.mat
      rm results/*.mat

      # fit permutation distribution
      echo matlab -nodesktop -nosplash -r "fitPerm;exit()"
      matlab -nodesktop -nosplash -r "fitPerm;exit()"

      # combine coherences across subjects
      cd ../stats
      echo matlab -nodesktop -nosplash -r "xsub_out_all_stats_allatl;exit()"
      matlab -nodesktop -nosplash -r "xsub_out_stats_allatl;exit()"

      # return to starting dir
      echo "cd $currDir"
      cd $currDir
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

