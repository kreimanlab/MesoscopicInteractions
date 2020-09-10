#!/usr/bin/python3

# Jiarui Wang :: jwang04@g.harvard.edu
# January 22, 2017

import sys
import os

# Color class
class col:
    H = '\033[95m' # Header
    B = '\033[94m' # Blue
    G = '\033[92m' # Green
    W = '\033[93m' # Warning
    E = '\033[91m' # Error
    e = '\033[0m' # ENDC

def getH5Files(pDir):
    h5fL = next(os.walk(pDir))[2]
    if ('.DS_Store' in h5fL):
        h5fL.remove('.DS_Store')
    H5files = []
    for h5fl in h5fL:
        if (h5fl.endswith('.hdf5')):
            H5files.append(h5fl)
    return H5files

def getConfig():
    # Open config file
    try:
        f = open('config','r')
        fl = f.readline()
        pdir = fl.split()[1]
        print('[*] Found H5DIR in config file: ' + pdir)

        # Remove trailing slash
        if (pdir.endswith('/')):
            pdir = pdir[:-1]

        return pdir
    except FileNotFoundError:
        print('[!] Error: config file not found.')
        exit()

#
# Prints the welcome splash art
#
def printSplash():
    print('#'+70*'-'+'#')
    print('|   Artifact Removal'+51*' '+'|')
    print('#'+70*'-'+'#')

def main():
    if (len(sys.argv) == 2):
        printSplash()
        pdir = getConfig()
        pid = str(sys.argv[1])
        pDir = pdir+'/'+pid
        h5files = getH5Files(pDir)
        print('[$] Starting artifact removal in: '+pDir)
        for hfile in h5files:
            print('\t(*) Using file: ', hfile)
            hfile = hfile.replace(' ','\ ')
            shCmd = 'matlab -nodesktop -nosplash -r "markmain('+"'"+pDir+'/'+hfile+"'"+')"'
            print(shCmd)
            os.system(shCmd)
    else:
        print(col.E+'[!] '+col.e+'Usage: python3 main.py m0000*')
    print('[!] All Done.')

if __name__ == '__main__':
    main()
