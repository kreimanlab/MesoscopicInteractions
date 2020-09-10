#
#      ~= Peter the Data Handler =~
#
#                \|||||/
#               ( ~   ~ )
#              @( 0   0 )@
#               (   C   )
#                \ \_/ /
#                 |___| 
#
#   I look through .txt raw data files to make sure everything is consistent and
#   provide a data summary (summary.txt) containing the size of the data
#

# TODO: firstlen might not be a good measure of the standard line length. Need 
# to modify to mode of first n lines to rule out rare chance first line is funky

import os
import sys
import time

# Color class
class col:
    H = '\033[95m' # Header
    B = '\033[94m' # Blue
    G = '\033[92m' # Green
    W = '\033[93m' # Warning
    E = '\033[91m' # Error
    e = '\033[0m' # ENDC

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
 
    return False

# takes the list of the current sample and checks it for any AMPSAT or other
# data saturation markers and fills it with values from the previous sample
def convertNumeric(ListCurr, PrevList, fn, linenum):
    List = ListCurr[2:-1]
    for i in range(len(List)):
        if (not is_number(List[i])):
            if (PrevList == []):
                #print('[!] Warn: non-num "{0}" in: {1}, line {2}. Prev samp doesnt exist. Writing value to zero!'.format(List[i], fn, linenum), file=sys.stderr)
                ListCurr[i+2] = '0';
            else:
                #print('[!] Warn: non-num "{0}" in: {1}, line {2}. Set to prev samp.'.format(List[i], fn, linenum), file=sys.stderr)
                ListCurr[i+2] = PrevList[i+2]
    return ListCurr

# parallelize on this function
def convertCSV(Flist_0):
    fn = Flist_0 + '.txt'
    foutlist = Flist_0 + '.csv'

    # Read in files one at a time and convert to csv 
    #print('[*] Converting txt to csv')
    #i = 0
    #Fs_all = []
    #for fn in Flist:
    # single-threaded version:
    #print('\twriting: ', Foutlist[i], '...')
    print('\twriting: ', foutlist, '...')
    fi = open(fn, 'r')
    fo = open(foutlist, 'w')
    s = 0
    firstlen = 0
    linenum = 1
        
    # sampling frequency check init
    tstamp = ''
    k = 0
    fs = 1
    Fs = []
    Lprev = []
    flens = []
        
    for line in fi:
        L = line.split()

        # Check for breaks in data
        if (L[0] == '---'):
            linenum = linenum - 1
            print('[!] Break in data in: {0}, line {1}. Skipping this line in csv conversion.'.format(fn,linenum),file=sys.stderr)
        else:

            # Check sampling frequency
            if (linenum == 1):
                # sample freq check
                tstamp = L[1]
            else:
                # sample freq check
                if (not tstamp == L[1]):
                    Fs.append(fs)
                    fs = 1
                    tstamp = L[1]
                else:
                    fs = fs + 1
            

            if (s < 10): # number of lines to consider for calculating expected line width
                #firstlen = len(L)
                flens.append(len(L))
                firstlen = max(set(flens), key=flens.count) # hacked to emulate mode(flens)
                s = s + 1

            # Check to make sure all lines are same length
            try:
                assert(len(L) == firstlen)
            except (AssertionError):
                print('[!] Warn: not all lines in data are same length.')
                print('[!] Warning in peter.py: not all lines in data file {0} are the same length (line {1})'.format(fn,linenum), file=sys.stderr)
                # Try to rescue extra or missing fields
                if (len(L) < firstlen):
                    print(col.H+'\tLine contains fewer fields than normal. Filling empty fields with zeros..'+col.e)
                    for ifill in range(firstlen-len(L)):
                        L.append('0')
                elif (len(L) > firstlen):
                    print(col.H+'\tLine contains more fields than normal. Trying to find trouble field..'+col.e)
                    try:
                        eval(L[2])
                        print(col.B + '\t(*) Elapsed time found where expected, trimming from the end..' + col.e)
                        for idel in range(len(L)-firstlen):
                            print('\t\t'+L.pop(len(L)-1))
                    except:
                        print(col.B + '\t(*) found: ' + L[2] + ' where elapsed time was expected. Deleting item from line, then trimming from end..' + col.e)
                        L.pop(2)
                        for idel in range(len(L)-firstlen):
                            print('\t\t'+L.pop(len(L)-1))
                    
                #exit()
        
            # convert list to all numeric values
            Lconv = convertNumeric(L, Lprev, fn, linenum)

            # write final list to file
            #print(str(Lconv))
            print(','.join(Lconv[2:-1]),file=fo)
            Lprev = Lconv
            
            # keep track of the line number
            
        linenum = linenum + 1
        
    #print(fn, Fs) # diagnostic only
    # save first Fs calculation
    try:
        Fs1 = Fs[0]
    except IndexError:
        Fs1 = fs

    Fs_sav = Fs
    Fs = set(Fs[1:len(Fs)]) # remove the very first Fs calculation
    try:
        FsRange = max(Fs) - min(Fs)
        if (FsRange > 2):
            print(col.E + '[!] Warn: file {0} contains variance in sampling frequency greater than 2 Hz'.format(fn) + col.e)
            print(col.E + 'Fs: {0}'.format(Fs) + col.e)

    except:
        #print('[!] Error, removed Fs:')
        #print(Fs)
        print('[!] Warn: sampling rate was not determined for file {0}'.format(fn))
        print(Fs_sav)
        #exit()
    
        #if (len(Fs) > 1):
        #print('[!] Error: file {0} contains inconsistent sampling frequency'.format(fn),file=sys.stderr)

    # single threaded version
    #Fs_all.append(Fs.pop())

    # multi threaded version 
    try:
        Fs_all = Fs.pop()
    except KeyError:
        Fs_all = Fs1

    #i = i + 1
    fi.close()
    fo.close()

    fsamp = open('peter_samples','a')
    print(Flist_0,linenum - 1,file=fsamp)
    fsamp.close()

    return Fs_all
 

def main():
    # handle input args
    if (len(sys.argv) == 2):
        efile = sys.argv[1]
        print('[*] Using exported file: '+efile)
    else:
        print('[!] Usage: python3 peter.py _Export-XXX.txt')
    
    # Get sampling rate directly from exported file
    ei = 1
    for line in open(efile,'r'):
        lineS = line.split()
        if (len(lineS) >= 5):
            line1 = lineS[1]
            line2 = lineS[2]
            if ((line1 == 'Sampling') and (line2=='Rate')):
                #FS = round(eval(lineS[3]))
                FS = eval(lineS[3])
                FS_units = lineS[4]
                print('[@] Sampling Rate: ',str(FS),' ',FS_units)
                break

        if (ei == 12):
            print('[!] Error: Sampling Rate not found in exported file.')
            exit()

        ei = ei + 1

    # configuration
    summaryfname = 'peter_summary'
    flistfname = 'peter_filelist'

    # welcome message
    wd = os.getcwd()
    dt = time.strftime("%Y-%m-%d, %H:%M:%S %p")
    print('o~---------------------------------------------------------------~o')
    print('    Hi! I am Peter the Data Handler.\n')
    print('    The current time is: ' + dt)
    print('    I am in: ' + wd)
    print('o~---------------------------------------------------------------~o')

    # get names of all .txt files in current directory
    print('[*] Reading in .txt files in current directory...')
    names = []
    for f in os.listdir('.'):
        if f.endswith('.txt'):
            names.append(f)

    # grab the prefix, series information, and suffix
    prefix = []
    series = []
    suffix = []
    for f in names:
        ff = f.split('_')
        try:
            prefix.append('_'.join(ff[0:-2]))
            series.append(int(ff[-2]))
            suffix.append(ff[-1])
        except (ValueError,IndexError):
            print('[!] Error: names of .txt files in this dir do not conform to the right format. Check that there are no .txt files in this directory that are not supposed to be there. Check that all .txt files are named: <patient identifier>_<patient identifier>_<data series>_<number of samples>.txt')
            exit()

    # check that all prefixes are identical
    print('[*] Checking that all file prefixes are identical...')
    for i in range(len(prefix)):
        try:
            assert(prefix[i] == prefix[0])
        except (AssertionError):
            print('[!] Error: not all file prefixes are identical. Check that all the data in this directory are from one patient.')
            exit()

    # check that all suffixes are identical
    print('[*] Checking that all file suffixes are identical...')
    for i in range(len(suffix)):
        try:
            assert(suffix[i] == suffix[0])
        except (AssertionError):
            print('[!] Error: not all file suffixes are identical. Check that all the data in this directory are from one patient.')
            exit()

    series = sorted(series)

    # check that there are no breaks in the series
    #print('[*] Checking that there are no breaks in the file series...')
    #SERIES_START = series[0]
    #SERIES_END = series[-1]
    #for i in range(SERIES_END - SERIES_START + 1):
    #    try:
    #        assert(series[i] == SERIES_START+i)
    #    except (AssertionError,IndexError):
    #        print('[!] Error: the series of the .txt files in this directory do not form a continuous series. Check to make sure you are not missing any files.')
    #        exit()

    # open up first file and extract useful information
    print('[*] Extracting file details... ')

    # check for sampling frequency
    tstamp = ''
    k = 0
    fs = 1
    Fs = []

    infile = open(prefix[0]+'_'+str(series[0])+'_'+suffix[0],'r')
    i = 1
    N_CHAN = 0;
    N_SAMPLES = 0;
    for line in infile:
        L = line.split()
        if (L[0] == '---'):
            # break in data
            break

        if (i == 1):
            # sample freq check on first line
            tstamp = L[1]

            for j in range(len(L)):
                try:
                    float(L[j]) # Only count samples that are numeric
                    N_CHAN = N_CHAN + 1
                except ValueError:
                    if ((L[j] == 'AMPSAT') or (L[j] == 'SHORT')):
                        # Also count AMPSAT, SHORT events as channels
                        N_CHAN = N_CHAN + 1
                    else:
                        pass
        else:
            # sample freq check
            if (not tstamp == L[1]):
                Fs.append(fs)
                fs = 1
                tstamp = L[1]
            else:
                fs = fs + 1

        if (len(L) != 0):
            N_SAMPLES = N_SAMPLES + 1;
        i = i + 1
    
    Fs = set(Fs[1:len(Fs)]) # remove the very first and last Fs calculation
    if (len(Fs) > 1):
        print('[!] Warn: first file contains inconsistent sampling frequency',file=sys.stderr)

    print('\tN_SAMPLES: {0}, N_CHAN: {1}'.format(N_SAMPLES,N_CHAN))
    print('\tFs (Hz): ', Fs)
    infile.close()

    # saving directory information
    print('[*] Saving file information to: '+summaryfname+', '+flistfname)
    
    summaryf = open(summaryfname,'w')
    flistf = open(flistfname,'w')
    os.system('rm peter_samples')
    print('N_SAMPLES\t{0}'.format(N_SAMPLES),file=summaryf)
    print('N_CHAN\t{0}'.format(N_CHAN),file=summaryf)

    Flist = []
    Foutlist = []
    Flist_0 = []
    for i in range(len(prefix)):
        flistx = prefix[i]+'_'+str(series[i])+'_'+suffix[i]
        flistxOut = prefix[i]+'_'+str(series[i])+'_'+suffix[i].split('.txt')[0]+'.csv'
        flistxOutp = prefix[i]+'_'+str(series[i])+'_'+suffix[i].split('.txt')[0]
        Foutlist.append(flistxOut)
        Flist.append(flistx)
        Flist_0.append(flistxOutp)
        print(flistxOutp,file=flistf)


    # --- Begin parallel threading ---
    # Inputs: Flist, Foutlist

    # Single-threaded implementation
    #Fs_all = convertCSV(Flist, Foutlist)
   
    # Multi-threaded implementation
    
    #i = 0
    Fs_all = []
    
    #from joblib import Parallel, delayed
    import multiprocessing

    inputs = Flist_0
    num_cores = multiprocessing.cpu_count()
    print('[$] Parallel job started on {0} cores.'.format(num_cores))
    print('[*] Converting txt to csv')
    pool = multiprocessing.Pool(num_cores)
    Fs_all = pool.map(convertCSV, inputs)
    #print('poolint: ', poolint)
    #print('*poolint: ', *poolint)
    #Fs_all = zip(*poolint)
    #Fs_all = zip(*pool.map(convertCSV, inputs))
    #print('Fs_all: ', Fs_all)

    #Fs_all = Parallel(n_jobs=num_cores)(delayed(convertCSV)(i) for i in inputs)

    # --- End parallel threading ---
    # Outputs: Fs_all
   
    # Set Fs_all to pre-determined sampling frequency
    Fs_all = []
    Fs_all.append(FS)

    # re-sort to pre-parallelization order
    fsamp = open('peter_samples','r')
    Flist_samp = []
    for line in fsamp:
        Flist_samp.append((line.split()[0],int(line.split()[1])))
    Flist_samp.sort(key=lambda x: Flist_0.index(x[0]))
    fsamp.close()
    fsampout = open('peter_samples','w')
    for z in range(len(Flist_0)):
        print(Flist_samp[z][1], file=fsampout)
        #print(Flist_samp[z][0], Flist_samp[z][1], file=fsampout)
    fsampout.close()
   
    print('[*] saving Fs information to {0} ...'.format(summaryfname))
    if (len(set(Fs_all)) == 1):
        print('Fs\t{0}'.format(Fs_all[0]),file=summaryf)
    else:
        print('[!] Warning: sampling frequency is inconsistent in this data directory. Writing the most common sampling frequency to summary. See the following list of all sampling frequencies:',Fs_all,file=sys.stderr)
        print('Fs\t{0}'.format(max(set(Fs_all), key=Fs_all.count)),file=summaryf)

    print('[!] Done.')

if __name__ == '__main__':
    main()
