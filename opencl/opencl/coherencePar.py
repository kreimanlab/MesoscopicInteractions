import h5py
import numpy
from scipy import signal
from scipy import stats
import matplotlib.pyplot as plt
import multiprocessing as mp

                                                       
HG_S = 70                                                                    
HG_E = 200    

def coh_worker(a):
    v1 = a[0]
    v2 = a[1]
    fs = a[2]
    f, C2 = signal.coherence(v1,v2,fs=fs, nperseg=fs*2, detrend=False)
    return [f,C2]

def coherence(v, fs):

    # Frequency boundaries in Hz
    PLI_S = 55
    PLI_E = 65
    DEL_S = 0.5
    DEL_E = 3
    THE_S = 3
    THE_E = 8
    ALP_S = 8
    ALP_E = 12
    BET_S = 12
    BET_E = 25
    GAM_S = 25
    GAM_E = 100
    BRO_S = 0.5
    BRO_E = 125
    N_BANDS = 6

    n_samples,n_chan = v.shape
    n_comb = int(0.5*n_chan*(n_chan-1))
    r = numpy.zeros((N_BANDS,n_comb))

    # Parallel
    Fall = numpy.zeros((n_chan,n_chan))
    C2all = numpy.zeros((n_chan,n_chan))
    A = []
    for i in range(1, n_chan):
        for j in range(i+1, n_chan+1):
            a = []
            a.append(v[:,i-1])
            a.append(v[:,j-1])
            a.append(fs)
            A.append(a)
    pool = mp.Pool()
    results = [pool.apply_async(coh_worker, args=(A[x],)) for x in range(len(A))]
    output = [p.get() for p in results]

    c_comb = 0
    for i in range(1, n_chan):
        for j in range(i+1, n_chan+1):
            #f, C2 = signal.coherence(v[:,i-1],v[:,j-1],fs=fs, nperseg=fs*2, detrend=False)
            f = output[c_comb][0]
            C2 = output[c_comb][1]
            mask_del = ((f > DEL_S) & (f < DEL_E))
            mask_the = ((f >= THE_S) & (f < THE_E))
            mask_alp = ((f >= ALP_S) & (f < ALP_E))
            mask_bet = ((f >= BET_S) & (f < BET_E))
            #mask_gam = ((f >= GAM_S) & (f < GAM_E))
            mask_gam = ((f >= GAM_S) & (f < PLI_S)) | ((f >= PLI_E) & (f < GAM_E))
            mask_bro = ((f >= BRO_S) & (f < BRO_E))
            r[0,c_comb] = numpy.mean(C2[mask_del])
            r[1,c_comb] = numpy.mean(C2[mask_the])
            r[2,c_comb] = numpy.mean(C2[mask_alp])
            r[3,c_comb] = numpy.mean(C2[mask_bet])
            r[4,c_comb] = numpy.mean(C2[mask_gam])
            r[5,c_comb] = numpy.mean(C2[mask_bro])

            #print(r[:,c_comb])
            c_comb = c_comb + 1

    #plt.plot(f, C2)
    #plt.xlabel('Frequency (Hz)')
    #plt.ylabel('Magnitude Squared Coherence')
    #plt.show()

    return r


def coherence2(v, w, fs):

    # Frequency boundaries in Hz
    PLI_S = 55
    PLI_E = 65
    DEL_S = 0.5
    DEL_E = 3
    THE_S = 3
    THE_E = 8
    ALP_S = 8
    ALP_E = 12
    BET_S = 12
    BET_E = 25
    GAM_S = 25
    GAM_E = 100
    BRO_S = 0.5
    BRO_E = 125
    N_BANDS = 6

    n_samples,n_chan = v.shape
    n_comb = int(0.5*n_chan*(n_chan-1))
    r = numpy.zeros((N_BANDS,n_comb))


    # Parallel
    Fall = numpy.zeros((n_chan,n_chan))
    C2all = numpy.zeros((n_chan,n_chan))
    A = []
    for i in range(1, n_chan):
        for j in range(i+1, n_chan+1):
            a = []
            a.append(v[:,i-1])
            a.append(w[:,j-1])
            a.append(fs)
            A.append(a)
    pool = mp.Pool()
    results = [pool.apply_async(coh_worker, args=(A[x],)) for x in range(len(A))]
    output = [p.get() for p in results]


    c_comb = 0
    for i in range(1, n_chan):
        for j in range(i+1, n_chan+1):
            #f, C2 = signal.coherence(v[:,i-1],w[:,j-1],fs=fs, nperseg=fs*2, detrend=False)
            f = output[c_comb][0]
            C2 = output[c_comb][1]
            mask_del = ((f > DEL_S) & (f < DEL_E))
            mask_the = ((f >= THE_S) & (f < THE_E))
            mask_alp = ((f >= ALP_S) & (f < ALP_E))
            mask_bet = ((f >= BET_S) & (f < BET_E))
            #mask_gam = ((f >= GAM_S) & (f < GAM_E))
            mask_gam = ((f >= GAM_S) & (f < PLI_S)) | ((f >= PLI_E) & (f < GAM_E))
            mask_bro = ((f >= BRO_S) & (f < BRO_E))
            r[0,c_comb] = numpy.mean(C2[mask_del])
            r[1,c_comb] = numpy.mean(C2[mask_the])
            r[2,c_comb] = numpy.mean(C2[mask_alp])
            r[3,c_comb] = numpy.mean(C2[mask_bet])
            r[4,c_comb] = numpy.mean(C2[mask_gam])
            r[5,c_comb] = numpy.mean(C2[mask_bro])

            #print(r[:,c_comb])
            c_comb = c_comb + 1

    #plt.plot(f, C2)
    #plt.xlabel('Frequency (Hz)')
    #plt.ylabel('Magnitude Squared Coherence')
    #plt.show()

    return r

def main():
    h5fname = '/mnt/cuenap2/scripts/synth/out/sub8.h5'
    h5f = h5py.File(h5fname,'r')
    fs = int(h5f['/h5eeg/eeg'].attrs['rate'][0])
    n_chan = int(h5f['/h5eeg/eeg'].attrs['n_chan'][0])
    w = 60
    r_idx = 1000
    r_rows = int(w*fs)
    X = h5f['/h5eeg/eeg'][r_idx:(r_idx+r_rows),0:(0+n_chan)]
    #for i in range(n_chan):
    #    X[:,i] = stats.rankdata(X[:,i])

    r_idx = 71337
    Y = h5f['/h5eeg/eeg'][r_idx:(r_idx+r_rows),0:(0+n_chan)]
    #for i in range(n_chan):
    #    Y[:,i] = stats.rankdata(Y[:,i])


    r_coh = coherence(X, fs)
    print(r_coh)
    print(r_coh[0,:])

    r_coh2 = coherence2(X, Y, fs)
    print(r_coh2)

if __name__ == '__main__':
    main()
