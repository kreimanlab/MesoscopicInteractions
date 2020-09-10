#!/bin/bash

echo "_  _ ____ _  _ _  _"
echo "|\/| |___ |\ | |  |"
echo "|  | |___ | \| |__|"
echo ""

PS3=": "

select opt in "run config.sh" "run make.sh" "run clean.sh" "run preprocessing.sh" "run coherence.sh" "run figures.sh" quit; do

  case $opt in
    "run config.sh")
      echo "./config.sh"
      ./config.sh
      break
      ;;
    "run make.sh")
      echo "./make.sh"
      ./make.sh
      break
      ;;
    "run clean.sh")
      echo "./clean.sh"
      ./clean.sh
      break
      ;;
    "run preprocessing.sh")
      echo "./preprocessing.sh"
      ./preprocessing.sh
      break
      ;;
    "run coherence.sh")
      echo "./coherence.sh"
      ./coherence.sh
      break
      ;;
    "run figures.sh")
      echo "./figures.sh"
      ./figures.sh
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

