#!/usr/bin/python3

# Jerry Wang :: jwang04@g.harvard.edu
# January 20, 2017

import sys
import os

def distillFname(name):
    return name.split('_')[-1]

def main():
    # Input
    try:
        pid = str(sys.argv[1])
    except IndexError:
        print('[!] Usage: python3 clearNames.py mXXXXX')
        exit()
    files = next(os.walk(pid))[2]

    # filter file list for hdf5 files
    for fname in files:
        if (fname.endswith('.hdf5')):
            newname = pid + '_' + distillFname(fname);
            #os.system('mv '+pid+'/'+fname+' '+pid+'/'+newname)
            shcmd = 'mv '+pid+'/'+fname.replace(' ','\ ')+' '+pid+'/'+newname
            print(shcmd)
            os.system(shcmd)

if __name__ == '__main__':
    main()
