#!/bin/bash
echo source downloads/pyenv/bin/activate
source downloads/pyenv/bin/activate

# Check for matlab
echo -e "\e[36m[*] Testing MATLAB install..\e[39m"
RESULT_MAT=`whereis matlab | tail -n 1`
if [ -z "$RESULT_MAT" ]
then
    echo "MATLAB is NOT installed."
    echo -e "\e[31mPlease install MATLAB.\e[39m"
else
    echo -e "found: \e[32m$RESULT_MAT\e[39m"
    matlab -nodesktop -nosplash -r "vers=version;disp(vers);fprintf('\nmagic(9):\n');disp(magic(9));exit();"
    echo -e "\e[33m[*] Please check that matlab is activated and running.\e[39m"
fi

# Check for python packages
echo -e "\e[36m[*] Looking for python3 packages using pip3..\e[39m"
for packagename in scipy numpy matplotlib h5py reikna pyopencl
do
    RESULT_PY=`pip3 list --format=legacy | grep -F $packagename`
    #echo $RESULT_SCIPY
    if [ -z "$RESULT_PY" ]
    then
          echo -e "\e[31m$packagename is NOT installed.\e[39m"
    else
          echo -e "found: \e[32m$RESULT_PY\e[39m"
    fi
done

echo "[!] If you use another python package manager (e.g. conda, etc), make sure the above packages are installed."

echo -e "\e[36m[!] Config script done. Please inspect output and install missing packages.\e[39m"
