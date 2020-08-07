# Jiarui "Jerry" Wang :: jwang04@g.harvard.edu
# Kreiman Lab :: klab.tch.harvard.edu

import sys
import pathlib
import socket

def main():

    path = get_paths(sys.argv)

def get_paths(argv):

    path = []
    # check input arguments
    if (len(sys.argv) == 4):
        path.append(sys.argv)
        path.append('graph')
    elif (len(sys.argv) == 5):
        path.append(sys.argv)
        path.append('perm')
    else:
        # usage message
        print('Usage: fconnect.py h5file metric [n_perm] dev\n')
        print('\th5file\tPath of .h5 file containing multivariate timeseries')
        print('\tmetric\tCorrelation metric to use\n')
        print('\t\tp\tPearson Correlation Coefficient')
        print('\t\ts\tSpearman Correlation Coefficient')
        print('\t\tpP\tpartial p')
        print('\t\tsP\tpartial s')
        print('\t\tpc\tmagnitude-square coherence')
        print('\t\tsc\tspearman magnitude-Square coherence')
        print('\t\tpe\tp of the filtered envelope')
        print('\t\tse\ts of the filtered envelope\n')
        print('\tn_perm\t[optional] compute this number of permutations')
        print('\tdev\t OpenCL device number to use\n')
        print('Examples:\n\nCompute the full correlation matrix:\n')
        print('\tfconnect.py m00001.h5 p 0\n')
        print('Compute time-shifted correlations with 10000 permutations:\n')
        print('\tfconnect.py m00001.h5 p 10000 0\n')
        exit()

    host = socket.gethostname()
    if (('Leibniz' in host) or ('seldon' in host)):
        out_dir = './results'
        h5fname = './h5/' + argv[1] + '.h5'
        if (not pathlib.Path(h5fname).is_file()):
            h5fname = '/Volumes/cuenap_ssd/h5/' + argv[1] + '.h5'
    elif ('o2.rc.hms.harvard.edu' in host):
        out_dir = '/n/scratch2/jw324/opencl/results'
        h5fname = '/n/scratch2/jw324/data/h5/' + argv[1] + '.h5'
    else:
        out_dir = '.'
        h5fname = './h5/' + argv[1] + '.h5'
        if (not pathlib.Path(h5fname).is_file()):
            h5fname = '/mnt/cuenap2/scripts/synth/out/' + argv[1] + '.h5'
    path.append(out_dir)
    path.append(h5fname)
    return path

if __name__ == "__main__":
    main()
