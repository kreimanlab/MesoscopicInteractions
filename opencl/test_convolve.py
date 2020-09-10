import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.ndimage import gaussian_filter
from scipy.ndimage import gaussian_filter1d

sid = 'sub3'
h5fname = '/media/jerry/internal/data/h5_notch20/' + sid + '.h5'

h5f = h5py.File(h5fname,'r')
n_chan = int(h5f['/h5eeg/eeg'].attrs['n_chan'][0])
n_samples = int(h5f['/h5eeg/eeg'].attrs['n_samples'][0])
fs = h5f['/h5eeg/eeg'].attrs['rate'][0]
w = 10
r_rows = int(w * round(fs))
r_idx = 20000
n_chan = 1
X = h5f['/h5eeg/eeg'][r_idx:(r_idx+r_rows),0:(0+n_chan)]
t = np.linspace(0,w,num=r_rows)
print(X.shape)
Xf = gaussian_filter1d(X[:,0], sigma=int(fs/10))

# plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t,X)
ax.plot(t,Xf)
ax.set(xlabel='Time (s)', ylabel='IFP (uV)')
plt.show()
