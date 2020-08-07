import numpy
import os

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def main():
    metrics = ['s','p','scDe','scTh','scAl','scBe','scGa','scBr','pcDe','pcTh','pcAl','pcBe','pcGa','pcBr','sP','pP','sd','st','sa','sb','sg']
    ftypes = ['graph','perm','dists']
    resultsDir = '/mnt/cuenap2/data/results'
    for path,dirs,files in os.walk(resultsDir):
        break

    # Build patient list
    pat = []
    for f in files:
        pat.append(f.split('_')[0])
    pat = numpy.unique(pat)

    # Print patient list
    numer = 0
    denom = 0
    cc = 0
    for p in pat:
        cc = cc + 1
        print('({1})\t{0}'.format(p,cc))
        fp = []
        for f in files:
            # Make patient specific file list
            if (f.startswith(p)):
                fp.append(f)

        # Print each metric whether a file type exists or not
        for m in metrics:
            reportA = m + ':\t'
            for t in ftypes:
                report = t.upper()
                #print(fp)
                #print(' '.join(fp))
                isFile = False
                for f2 in fp:
                    if ((t in f2) and (('-'+m) in f2)):
                        isFile = True
                    
                #if (t in ' '.join(fp)) and (m in ' '.join(fp)):
                if (isFile):
                    report = bcolors.OKGREEN + report + bcolors.ENDC
                    numer = numer + 1
                else:
                    report = bcolors.FAIL + report + bcolors.ENDC
                report = report + ' '
                reportA = reportA + report
                denom = denom + 1
            print('\t\t{0}'.format(reportA))

    print('[*] Total number of patients: {0}'.format(len(pat)))
    print('[*] Percentage completed: {0:0.2f}%'.format(100*numer/denom))

main()
