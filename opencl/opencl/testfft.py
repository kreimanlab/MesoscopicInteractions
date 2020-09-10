
import pyfftw
import numpy
from timeit import time
from scipy.signal import hilbert

m = 250*60
n = 128
a = numpy.random.normal(0,1,(n,m))

t = time.time()
a2 = pyfftw.interfaces.numpy_fft.fft(a, axis=1)
#b = hilbert(a,axis=0)
#c = abs(b)
print(time.time()-t)

t = time.time()
#a2 = pyfftw.interfaces.numpy_fft.fft(a, axis=0)
a2 = numpy.fft.fft(a, axis=1)
#b = hilbert(a,axis=0)
#c = abs(b)
print(time.time()-t)

