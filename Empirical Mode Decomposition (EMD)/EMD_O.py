import numpy as np
import scipy.io as sio
from numpy import sqrt, square
from scipy.signal import argrelmax, argrelmin
import matplotlib.pyplot as plt
from scipy import interpolate

from PyEMD import EMD
import numpy as np

fs = 500
f = 5
sample = 500
x = np.arange(sample) #x is array number of samples from 0 to maxnumber of samples

s = np.sin(2*np.pi*f*x/fs) + np.sin(10*np.pi*f*x/fs)

emd = EMD()
IMFs = emd(s)

print(len(IMFs))


f, axarr = plt.subplots(7, sharex=True)
axarr[0].plot(s, 'b')
axarr[1].plot(IMFs[0], 'r')
axarr[2].plot(IMFs[1], 'r')
axarr[3].plot(IMFs[2], 'r')
axarr[4].plot(IMFs[3], 'r')
axarr[5].plot(IMFs[4], 'r')
axarr[6].plot(IMFs[0]+IMFs[1]+
              IMFs[2]+IMFs[3]+
              IMFs[4], 'r')

plt.show()

