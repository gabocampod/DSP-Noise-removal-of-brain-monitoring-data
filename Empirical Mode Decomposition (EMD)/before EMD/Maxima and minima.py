import numpy as np
from scipy.signal import argrelmax
from scipy.signal import argrelmin
import matplotlib.pyplot as plt
from scipy import interpolate

fs = 500
f = 5
sample = 500
x = np.arange(sample)

s = np.sin(2*np.pi*f*x/fs)   + np.cos(10*np.pi*f*x/fs)

maxima_arr = argrelmax(s)

minima_arr = argrelmin(s)

print(maxima_arr)
print(minima_arr)

print(s[argrelmax(s)[0]])
print(s[argrelmin(s)[0]])


plt.figure()
plt.plot(x, s, 'r', argrelmax(s)[0], s[argrelmax(s)[0]], '--b', argrelmin(s)[0], s[argrelmin(s)[0]], '--g')
plt.show()