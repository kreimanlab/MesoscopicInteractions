import os
import re

def tryint(s):
    try:
        return int(s)
    except:
        return s
             
def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)
    return l

def cleanOffs(f):
    if f.endswith('.txt'):
        print('[$] Converting '+f)
        if (len(f.split('_')) == 5):
            fn = f.split('_')
            fn.pop(2)
            fn = '_'.join(fn)
            #print('\t\trenaming '+f+' to '+fn)
            os.system('mv '+f+' '+fn)
            f = fn
        else:
            pass
            #print('\t\tText filenames conform to format. Not renaming.')
        # Remove Offs if found
        infi = open(f,'r')
        newf = open('new'+f,'w')
        for line in infi:
            if (line.split()[2] == 'OFF'):
                nl = line.split()
                nl.pop(2)
                newl = nl[0]+' '+nl[1]+'\t'+'\t'.join(nl[2:])
                print(newl,file=newf)
            else:
                print(line.split('\n')[0],file=newf)
        infi.close()
        newf.close()
        os.system('mv new'+f+' '+f)


def main():
    os.system('rm new*.txt')
    Files = next(os.walk('.'))[2]
    
    # --- Begin parallel threading ---
    # Inputs: Files

    # Single-threaded implementation
    #for f in sort_nicely(Files):
    #    cleanOffs(f)
    
    # Multi-threaded implementation
    import multiprocessing

    inputs = Files
    num_cores = multiprocessing.cpu_count()
    print('[$] Parallel job started on {0} cores.'.format(num_cores))
    print('[*] Removing OFF from data...')
    pool = multiprocessing.Pool(num_cores)
    Fs_all = pool.map(cleanOffs, inputs)

    # --- End parallel threading ---
    # Outputs: null


    print('[*] Cleaning up extra new files..')
    os.system('rm new*.txt')
    print('[!] All done.')


if __name__ == '__main__':
    main()
