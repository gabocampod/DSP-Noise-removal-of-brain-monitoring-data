import numpy as np
from numpy import sqrt, square
from scipy.signal import argrelmax
from scipy.signal import argrelmin
import matplotlib.pyplot as plt
from scipy import interpolate


Fs = 8000
f = 5
sample = 800
x = np.arange(sample)

s2 =np.random.rand(80)
s = np.sin(2* np.pi * f *  x / Fs) + np.cos(100 * np.pi * f * x / Fs) + np.cos(60 * np.pi * f * x / Fs)


j = 0

component = np.zeros(shape=(x.size,10))

data_to_emd = s



while argrelmax(data_to_emd)[0].size > 3 and max(abs(data_to_emd)) > 1e-10:


    ###MAXIMA AND MINIMA CALCULATIONS
    maxima_arr = argrelmax(data_to_emd)
    minima_arr = argrelmin(data_to_emd)

    tck_max = interpolate.splrep(maxima_arr[0], data_to_emd[maxima_arr[0]], s=0)
    tck_min = interpolate.splrep(minima_arr[0], data_to_emd[minima_arr[0]], s=0)

    max_interpolated = interpolate.splev(x, tck_max, der=0)
    min_interpolated = interpolate.splev(x, tck_min, der=0)


    ##MEAN

    mean = (max_interpolated + min_interpolated) / 2

    # prints mean as before but also mean
    plt.figure()
    plt.plot(x, data_to_emd, 'r', x, max_interpolated, '--b', x, min_interpolated, '--g', x, mean, 'k')
    plt.show()

    ### ELEMENT H1 CALCULATION

    h1 = data_to_emd - mean

    h_old = h1



    sd = 1.0

    while sd > 0.2:

        ###CUBIC SPLINE CALCULATIONS

        maxima_arr = argrelmax(h_old)
        minima_arr = argrelmin(h_old)

        tck_max = interpolate.splrep(maxima_arr[0], h_old[maxima_arr[0]], s=0)
        tck_min = interpolate.splrep(minima_arr[0], h_old[minima_arr[0]], s=0)

        max_interpolated = interpolate.splev(x, tck_max, der=0)
        min_interpolated = interpolate.splev(x, tck_min, der=0)


        ##MEAN CALCULTIONS
        mean = (max_interpolated + min_interpolated)/2

        #prints mean as before but also mean
        #plt.figure()
        #plt.plot(x, s, 'r', x, max_interpolated, '--b', x, min_interpolated, '--g', x, mean, 'k')
        #plt.show()

        ### ELEMENT CALCULATION

        h_new = h_old - mean

        #prints mean as before but also mean
        #plt.figure()
        #plt.plot(x, h_old, 'r', x, max_interpolated, '--b', x, min_interpolated, '--g', x, h_new, 'k',x, mean, '--y')
        #plt.show()

        sd = 0 #re-star SD on each loop

        for i in range (0,x.size,1):

            sd = sd + (pow(abs(h_old[i] - h_new[i]), 2))/(pow(h_old[i], 2))

            print(sd)

        h_old = h_new

    component[:,j] = h_old
    data_to_emd = data_to_emd - component[:,j]

    ax1 = plt.subplot(211)
    ax1.plot(data_to_emd,'b')

    ax2 = plt.subplot(212, sharex=ax1)
    ax2.plot(component[:,j],'r')

    plt.show()

    print(argrelmax(data_to_emd)[0].size)

    j = j+1

#end outter while loop


residual = data_to_emd[:]

sum_com = 0

for m in range(j):

    sum_com = sum_com + component [:,m]

initial_reconstructed =  sum_com + residual

comparisson = s - initial_reconstructed

rms = sqrt(np.mean(square(s-initial_reconstructed)))

f, axarr = plt.subplots(j+4, sharex=True)
axarr[0].plot(s,'b')

for k in range(1,j+1,1):
    axarr[k].plot(component[:,k-1],'r')


axarr[j+1].plot(residual,'b')
axarr[j+2].plot(initial_reconstructed,'b')
axarr[j+3].plot(comparisson,'b')
plt.show()

print(rms)
