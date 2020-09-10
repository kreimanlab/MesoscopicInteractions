# Last edited: January 7, 2017
# Jiarui Wang :: jwang04@g.harvard.edu

import os
import sys

def main():
    # Handle user input
    if (len(sys.argv) == 2):
        infilename = sys.argv[1]
    else:
        print('[!] Usage: python3 chunkraw.py [input_file]')
        exit()
    #print('[!JESUS] '+infilename)
    #exit()

    SKIP_LINES = 15
    LINES_PER_FILE = 100000
    i = 1
    j = 1
    currDir = os.getcwd().split('/')[-1]
    
    ## Override naming convention
    currDir = 'AA_00000000'

    #infile = open('_Export-TC_4ce6814b_Sleep_1.txt','r')
    infile = open(infilename,'r')
    for k in range(SKIP_LINES):
        infile.readline()

    outfname = currDir+'_'+str(j)+'_'+str(LINES_PER_FILE)+'.txt'
    outfile = open(outfname,'w')
    print('[!] Opened {0} for writing.'.format(outfname))
    for line in infile:
        # print only lines after SKIP_LINES
        #if (i > SKIP_LINES):
        print(line.split('\n')[0],file=outfile)

        if (i%LINES_PER_FILE == 0):
            j = j + 1
            outfile.close()
            outfname = currDir+'_'+str(j)+'_'+str(LINES_PER_FILE)+'.txt'
            outfile = open(outfname,'w')
            print('[!] Opened {0} for writing.'.format(outfname))

        #if (i == SKIP_LINES+250000):
        #    break

        i = i + 1
    infile.close()

main()
