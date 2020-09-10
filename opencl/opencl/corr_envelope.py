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

def corr_envelope_gpu(Xbip, thr):
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
    # --------------------------------------------------------------------------
    # Frequency boundaries in Hz
    #
    #   Any edits here requires also editing the frequency segmentation code
    #   immediately following coherence estimate (see below)
    #
    DEL_S = 0.5     # Delta wave
    DEL_E = 3
    THE_S = 3       # Theta wave
    THE_E = 8
    ALP_S = 8       # Alpha wave
    ALP_E = 12
    BET_S = 12      # Beta wave
    BET_E = 25
    GAM_S = 25      # Gamma wave
    GAM_E = 100
    BRO_S = 0.5     # Broadband
    BRO_E = 125
    N_BANDS = 6     # Number of frequency bands (not including power line)
    # --------------------------------------------------------------------------


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

   


def main():
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


if __name__ == '__main__':
    main()
