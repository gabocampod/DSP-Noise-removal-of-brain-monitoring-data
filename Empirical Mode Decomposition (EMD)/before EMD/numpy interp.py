import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.signal import argrelmax
from scipy.signal import argrelmin
from numba import njit
import cProfile

@njit
def _my(maxima_arr,minima_arr):

    str = 'cubic'
    str2 = 'extrapolate'

    max_interpolated = interpolate.interp1d(maxima_arr[0], s[maxima_arr[0]], kind=str, fill_value=str2)

    max_interpolated(x)


Fs = 500
f = 5
sample = 5000000

x = np.arange(sample)

s = np.sin(2 * np.pi * f * x / Fs) #+ np.cos(10 * np.pi * f * x / Fs) + np.sin(20 * np.pi * f * x / Fs)

maxima_arr = argrelmax(s)

minima_arr = argrelmin(s)

_my(maxima_arr,minima_arr)


