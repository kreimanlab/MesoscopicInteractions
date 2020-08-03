#!/usr/bin/python3

# Jiarui Wang :: jwang04@g.harvard.edu
# December 20, 2016

import os
import sys
from converter import *

#
# Reads config file and makes one with default values if it does not exist
#
# Returns a 2-element vector of directory strings
#
def getConfig():
    # Open config file
    ConfigDir = []
    try:
        f = open('config','r')
        Ents = []
        Dirs = []
        for line in f:
            line = line.split('\n')[0]
            Ents.append(line.split('=')[0])
            Dirs.append(line.split('=')[1])
        print('[!] Found following entries in config file:')
        ConfigDir = Dirs[0:2]
        print('\t' + Ents[0] + ': ' + Dirs[0])
        print('\t' + Ents[1] + ': ' + Dirs[1])
        if (len(ConfigDir) != 2):
            print('[!] Error: number of entries in config file is not 2.')
    except FileNotFoundError:
        print('[!] Warning: config file not found, setting default values.')
        ConfigDir.append('../exported')
        ConfigDir.append('../finished')
        Ents = ['EXPORTED_DIR','FINISHED_DIR']
        Dirs = ConfigDir[0:2]
        print('\t' + Ents[0] + ': ' + Dirs[0])
        print('\t' + Ents[1] + ': ' + Dirs[1])
        # Write config file with defaults
        f = open('config','w')
        for i in range(2):
            print(Ents[i]+'='+Dirs[i],file=f)
    # Strip terminal slash
    for j in range(2):
        if Dirs[j].endswith('/'):
            Dirs[j] = Dirs[j][0:-1]
    # Return final value
    return ConfigDir

#
# Prints the welcome splash art
#
def printSplash():
    print('#'+70*'-'+'#')
    print('|   Text to H5eeg conversion'+43*' '+'|')
    print('#'+70*'-'+'#')

#
# Main program
#
def main():
    # Check for input args
    if (len(sys.argv) == 2):
        printSplash()
        Dirs = getConfig()
        p1 = converter(Dirs, sys.argv[1])
        p1.buildExportList()
        p1.chunkExported()
        print(col.G+'[!] All done.'+col.e)
        #p1.printContents()
    else:
        print(col.W+'[!] '+col.e+'Usage: python3 main.py m0000*')

if __name__ == '__main__':
    main()
