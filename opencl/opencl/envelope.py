# Jiarui "Jerry" Wang :: jwang04@g.harvard.edu
# Kreiman Lab :: klab.tch.harvard.edu

import scipy
import h5py
from scipy.signal import butter, lfilter, hilbert
import matplotlib.pyplot as plt
import numpy
import reikna
from reikna.cluda import dtypes, any_api
from reikna.fft import FFT
from reikna.core import Annotation, Type, Transformation, Parameter
from timeit import time

# Profiling results
# get api:  4.76837158203125e-07
# init device:  2.765655517578125e-05
# type cast:  0.0010349750518798828
# real2imag:  0.001692056655883789
# init:  0.0027701854705810547
# compile1:  0.04964756965637207
# cpy to device:  0.00016021728515625
# allocate result on dev:  4.887580871582031e-05
# compute:  0.005811929702758789
# cpy from device:  0.01568126678466797
# hilbert:  0.029632091522216797
# compile2:  0.04875349998474121
# ifft:  0.017459630966186523

def envelope(Xbip):
    return abs( hilbert(Xbip, axis=0) )

def envelope_gpu(Xbip, thr):
    #
    #   GPU implementation of estimating the envelope of multiple timeseries
    #   data based on scipy.signal.hilbert. Reports the real-valued envelope of
    #   a multivariate input
    #
    #   Inputs:
    #       Xbip        2-D matrix where each column is a variable and each row
    #                   is a time sample
    #
    #       thr         GPU thread object (reikna.cluda.api.Thread)
    #
    #                       import reikna
    #                       dev_n = 0 # GPU device index, this varies by machine
    #                       plat_n = 0 # platform index, this is usually 0
    #                       api = reikna.cluda.any_api()
    #                       plat = api.get_platforms()[plat_n]
    #                       p0_dev = plat.get_devices()
    #                       dev = p0_dev[dev_n] 
    #                       thr = api.Thread(dev)
    #
    #   Outputs:
    #       X           same dimensions as input Xbip, real-valued matrix
    #                   containing the envelope of Xbip
    #
    #   Testing:
    #       95 pct difference from scipy.signal.coherence: 0.000727467, n=15336
    #
    arr = Xbip.astype(numpy.float32);
    # ----------------------------------------------------- 0.0010349750518798828
    # Real to complex
    trf = get_complex_trf(arr)
    # Compile
    fft = FFT(trf.output, axes=(0,))
    fft.parameter.input.connect(trf, trf.output, new_input=trf.input)
    cfft = fft.compile(thr)
    # ----------------------------------------------------- 0.04964756965637207
    arr_dev = thr.to_device(arr)
    # ----------------------------------------------------- 0.00016021728515625
    res_dev = thr.array(arr.shape, numpy.complex64)
    # ----------------------------------------------------- 4.887580871582031e-05
    cfft(res_dev, arr_dev, 0)
    # ----------------------------------------------------- 0.005811929702758789
    result = res_dev.get()
    # ----------------------------------------------------- 0.01568126678466797
    # Hilbert transform coefficient
    h = numpy.zeros(Xbip.shape, numpy.complex64)
    h_rows,h_cols = Xbip.shape
    for i in range(h_cols):
        if h_rows % 2 == 0: # if number of time samples is even
            h[0,i] = 1
            h[int(h_rows/2),i] = 1
            h[1:int(h_rows/2),i] = 2
        else: # if number of time samples is odd
            h[0,i] = 1
            h[1:int((h_rows+1)/2),i] = 2

    arr = numpy.multiply(result,h)
    # ----------------------------------------------------* 0.029632091522216797
    # Compile
    ffti = FFT(arr, axes=(0,))
    cffti = ffti.compile(thr)
    # ----------------------------------------------------- 0.04875349998474121
    # Compute inverse
    arr_dev = thr.to_device(arr)
    res_dev = thr.array(arr.shape, numpy.complex64)
    cffti(arr_dev, arr_dev, 1)
    Yf_h = arr_dev.get()
    # ----------------------------------------------------- 0.017459630966186523
    return abs(Yf_h)



def envelope_gpu_test(Xbip, thr):
    #SHOW_INFO = False

    # -----------------------------------------------------
    #t0 = time.time()

    #api = any_api()
    #plat = api.get_platforms()[plat_n]
    #p0_dev = plat.get_devices()

    # ----------------------------------------------------- 4.76837158203125e-07
    #t0p1 = time.time()
    #print('get api: ', t0p1 - t0)

    #if (SHOW_INFO):
    #    print('\tPlatform:\t',plat.name)
    #    print('\tVendor: \t',plat.vendor)
    #    print('\tVersion:\t',plat.version)
    #    print('\tDevices:\t')
    #    dc = 0
    #    for dev in p0_dev:
    #        if (dc == dev_n):
    #            dstr = '*'
    #        else:
    #            dstr = ' '
    #        print('\t\t[{0}]'.format(dstr),dev.name)
    #        print('\t\t\tLocal:\t {0:0.0f} K'.format(dev.local_mem_size/1e3))
    #        print('\t\t\tGlobal:\t {0:0.0f} G'.format(dev.global_mem_size/1e9))
    #        print('\t\t\tCores:\t',dev.max_compute_units)
    #        print('\t\t\tClock:\t',dev.max_clock_frequency,'MHz')
    #        dc = dc + 1
    #dev = p0_dev[dev_n]

    #thr = api.Thread.create(True)
    #thr = api.Thread(dev)

    # ----------------------------------------------------- 2.765655517578125e-05
    #t0p2 = time.time()
    #print('init device: ', t0p2 - t0p1)
    
    #arr = numpy.random.normal(size=3000).astype(numpy.float32)
    arr = Xbip.astype(numpy.float32);

    # ----------------------------------------------------- 0.0010349750518798828
    #t0p3 = time.time()
    #print('type cast: ', t0p3 - t0p2)

    # destination imagination
    #arr = numpy.ones(arr.shape, numpy.float32)
    trf = get_complex_trf(arr)
    #trf = arr.astype(numpy.complex64)

    # ----------------------------------------------------- 0.001692056655883789
    #t0p4 = time.time()
    #print('real2imag: ', t0p4 - t0p3)

    # ----------------------------------------------------- 0.0027701854705810547
    #t1 = time.time()
    #print('init: ', t1 - t0)

    # Compile
    fft = FFT(trf.output, axes=(0,))
    fft.parameter.input.connect(trf, trf.output, new_input=trf.input)
    cfft = fft.compile(thr)

    # ----------------------------------------------------- 0.04964756965637207
    #t2 = time.time()
    #print('compile1: ', t2 - t1)
    #t3 = time.time()

    # Compute


    arr_dev = thr.to_device(arr)
    # ----------------------------------------------------- 0.00016021728515625
    #t3p1 = time.time()
    #print('cpy to device: ', t3p1 - t2)

    res_dev = thr.array(arr.shape, numpy.complex64)
    # ----------------------------------------------------- 4.887580871582031e-05
    #t3p2 = time.time()
    #print('allocate result on dev: ', t3p2 - t3p1)

    cfft(res_dev, arr_dev, 0)
    # ----------------------------------------------------- 0.005811929702758789
    #t3p3 = time.time()
    #print('compute: ', t3p3 - t3p2)

    result = res_dev.get()
    # ----------------------------------------------------- 0.01568126678466797
    #t3p4 = time.time()
    #print('cpy from device: ', t3p4 - t3p3)

    # Validate
    #reference = numpy.fft.fft(arr, axis=0)
    #print(result, result.shape)
    #print(reference, reference.shape)
    #verr = numpy.linalg.norm(result - reference) / numpy.linalg.norm(reference)
    #print(verr)

    # Hilbert transform coefficient
    h = numpy.zeros(Xbip.shape, numpy.complex64)
    h_rows,h_cols = Xbip.shape
    for i in range(h_cols):
        if h_rows % 2 == 0: # if number of time samples is even
            h[0,i] = 1
            h[int(h_rows/2),i] = 1
            h[1:int(h_rows/2),i] = 2
        else: # if number of time samples is odd
            h[0,i] = 1
            h[1:int((h_rows+1)/2),i] = 2

    #print('h_rows: ',h_rows,' h_cols: ',h_cols)
    #print('h:\n',h)
    
    arr = numpy.multiply(result,h)
    #result = numpy.multiply(result,h)
    #print('mult:\n',result)
    #print(type(result))
    #arr = result;
    
    # ----------------------------------------------------- 0.029632091522216797
    #t4 = time.time()
    #print('hilbert: ', t4 - t3)

    # Compile
    ffti = FFT(arr, axes=(0,))
    cffti = ffti.compile(thr)

    # ----------------------------------------------------- 0.04875349998474121
    #t5 = time.time()
    #print('compile2: ', t5 - t4)

    # Compute inverse
    arr_dev = thr.to_device(arr)
    res_dev = thr.array(arr.shape, numpy.complex64)
    cffti(arr_dev, arr_dev, 1)
    Yf_h = arr_dev.get()
    #Yf_h = numpy.fft.ifft(arr,axis=0)
    
    # ----------------------------------------------------- 0.017459630966186523
    #t6 = time.time()
    #print('ifft: ', t6 - t5)
    
    return abs(Yf_h)


def butter_bandpass(lowcut, highcut, fs, order):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def bandpassf(data, a, b, fs):
    #b, a = butter_bandpass(lowcut, highcut, fs, o)
    y = lfilter(b, a, data)
    return y

def get_complex_trf(arr):
    try:
        complex_dtype = dtypes.complex_for(arr.dtype)
    except KeyError:
        complex_dtype = 'complex64'
    return Transformation(
        [Parameter('output', Annotation(Type(complex_dtype, arr.shape), 'o')),
        Parameter('input', Annotation(arr, 'i'))],
        """
            ${output.store_same}(
            COMPLEX_CTR(${output.ctype})(
                ${input.load_same},
                0));
        """)


def hilbert_custom(Xf):
    return hilbert(Xf,axis=0)

def main():
    #h5fname = '/mnt/cuenap2/scripts/synth/out/sub8.h5'
    h5fname = './h5/sub11.h5'
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

    print(X.shape)
    DEL_S = 0.5
    DEL_E = 1.5
    b,a = butter_bandpass(DEL_S,DEL_E,fs,3)
    #for i in range(n_chan):
    c = 8
    Xf = bandpassf(numpy.transpose(X),a,b,fs)
    Yf = bandpassf(numpy.transpose(Y),a,b,fs)
    #Xf = Xf[c,:]
    Xf = numpy.transpose(Xf)
    Yf = numpy.transpose(Yf)


    # Hilbert transform

    # Device number
    dev_n = 0
    # Platform number
    plat_n = 0

    api = any_api()
    plat = api.get_platforms()[plat_n]
    print('\tPlatform:\t',plat.name)
    print('\tVendor: \t',plat.vendor)
    print('\tVersion:\t',plat.version)
    print('\tDevices:\t')
    p0_dev = plat.get_devices()
    dc = 0
    for dev in p0_dev:
        if (dc == dev_n):
            dstr = '*'
        else:
            dstr = ' '
        print('\t\t[{0}]'.format(dstr),dev.name)
        print('\t\t\tLocal:\t {0:0.0f} K'.format(dev.local_mem_size/1e3))
        print('\t\t\tGlobal:\t {0:0.0f} G'.format(dev.global_mem_size/1e9))
        print('\t\t\tCores:\t',dev.max_compute_units)
        print('\t\t\tClock:\t',dev.max_clock_frequency,'MHz')
        dc = dc + 1
    dev = p0_dev[dev_n]

    #thr = api.Thread.create(True)
    thr = api.Thread(dev)

    #arr = numpy.random.normal(size=3000).astype(numpy.float32)
    arr = Yf.astype(numpy.float32);

    # destination imagination
    #arr = numpy.ones(arr.shape, numpy.float32)
    trf = get_complex_trf(arr)
    #trf = arr.astype(numpy.complex64)

    # Compile
    fft = FFT(trf.output, axes=(0,))
    fft.parameter.input.connect(trf, trf.output, new_input=trf.input)
    cfft = fft.compile(thr)

    # Compute
    arr_dev = thr.to_device(arr)
    res_dev = thr.array(arr.shape, numpy.complex64)
    cfft(res_dev, arr_dev, 0)
    result = res_dev.get()

    # Validate
    reference = numpy.fft.fft(arr, axis=0)
    print(result, result.shape)
    print(reference, reference.shape)
    verr = numpy.linalg.norm(result - reference) / numpy.linalg.norm(reference)
    print(verr)

    # Hilbert transform coefficient
    h = numpy.zeros(Xf.shape, numpy.complex64)
    h_rows,h_cols = Xf.shape
    for i in range(h_cols):
        if h_rows % 2 == 0: # if number of time samples is even
            h[0,i] = 1
            h[int(h_rows/2),i] = 1
            h[1:int(h_rows/2),i] = 2
        else: # if number of time samples is odd
            h[0,i] = 1
            h[1:int((h_rows+1)/2),i] = 2

    print('h_rows: ',h_rows,' h_cols: ',h_cols)
    print('h:\n',h)
    
    result = numpy.multiply(result,h)
    print('mult:\n',result)
    print(type(result))
    arr = result;
    

    # Compile
    ffti = FFT(arr, axes=(0,))
    #fft.parameter.input.connect(trf, trf.output, new_input=trf.input)
    cffti = ffti.compile(thr)

    # Compute inverse
    arr_dev = thr.to_device(arr)
    res_dev = thr.array(arr.shape, numpy.complex64)
    cffti(arr_dev, arr_dev, 1)
    Yf_h = arr_dev.get()
    #Yf_h = numpy.fft.ifft(arr,axis=0)
    
    XenvNew = abs(Yf_h)


    # reference hilbert
    Xenv = abs(hilbert(Yf,axis=0))
    Xenva = abs(hilbert(Xf,axis=0))

    # fast hilbert
    Xenvb = envelope_gpu(Yf,thr)

    print('New hilbert: ',XenvNew[0:4,0:10],XenvNew.shape)
    print('Fast hilbert: ',Xenvb[0:4,0:10],Xenvb.shape)
    print('scipy hilbert: ',Xenv[0:4,0:10],Xenv.shape)
    Xenv = XenvNew

    # Calculate corr matrix
    C = numpy.cov(X,rowvar=False)
    Cd = numpy.sqrt(numpy.diagonal(C))

    # Test
    d = numpy.abs(Xenv - Xenvb).ravel()
    d = sorted(d)
    print('95 pct difference: ', d[int(0.95*len(d))])
    print('max difference: ', max(d))
    print('n=', len(d))

    # Plot
    Yf = Yf[:,c]
    Xf = Xf[:,c]
    Yenv = Xenv[:,c]
    Xenva = Xenva[:,c]
    T = numpy.linspace(0,len(Xf)/fs,len(Xf))
    print(Xf,len(Xf))
    print(Yenv,len(Yenv))
    print(Xenva,len(Xenva))
    print('Corr coef: ',numpy.corrcoef(scipy.stats.rankdata(Xf),scipy.stats.rankdata(Yf),rowvar=False))
    plt.plot(T,Xf,'black',T,Yf,'gray',T,Yenv,'orange',T,Xenva,'g',linewidth=1.0)
    plt.xlabel('Time (sec)')
    plt.ylabel('Voltage (uV)')
    plt.show()


    return


if __name__ == '__main__':
    main()
