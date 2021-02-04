import numpy as np
from scipy.signal import argrelmax
from scipy.signal import argrelmin
import matplotlib.pyplot as plt
from scipy import interpolate

def find_all_extrema(data_emd, x, max_array_positions, min_array_positions):

    #(signal,x sample values, max_arr[0],min_arr[0])

    nbsym = 2
    MIN_END, MAX_END = len(min_array_positions), len(max_array_positions)

    #LEFT BOUNDARY

    if  max_array_positions[0] < min_array_positions[0]:         #if maxima occurs before minima
        if data_emd[0] > data_emd[min_array_positions[0]]:                            #if data_emd[0] is even larger than the first maxima

            lmax = max_array_positions[1:min(MAX_END, nbsym + 1)][::-1] #[::-1] np.flipr but for 1D array
            lmin = min_array_positions[0:min(MIN_END, nbsym)][::-1]
            lsym = max_array_positions[0]
        else:                                                                         # if data_emd[0] is less than first maxima
            lmax = max_array_positions[0:min(MAX_END, nbsym)][::-1]
            lmin = np.append(min_array_positions[1:min(MIN_END, nbsym-1)][::-1] ,0)
            lsym = 0
    else:                                                       #if minima occurs before maxima
        if data_emd[0] < data_emd[max_array_positions[0]]:                            #if data_emd[0] is even less than first minima
            lmax = max_array_positions[0:min(MAX_END, nbsym)][::-1]
            lmin = min_array_positions[1:min(MIN_END, nbsym + 1)][::-1]
            lsym = min_array_positions[0]
        else:                                                                         #if data_emd[0] is more than first minima
            lmax = np.append(max_array_positions[0:min(MAX_END, nbsym - 1)][::-1], 0)
            lmin = min_array_positions[0:min(MIN_END, nbsym)][::-1]
            lsym = 0

    #RIGHT BOUNDARY

    if max_array_positions[-1] < min_array_positions[-1]:        #If last extrema is a minima
        if data_emd[-1] < data_emd[max_array_positions[-1]]:                          #If last value of data is less than last maxima
            rmax = max_array_positions[max(MAX_END - nbsym, 0):][::-1]
            rmin = min_array_positions[max(MIN_END - nbsym - 1, 0):-1][::-1]
            rsym =min_array_positions[-1]
        else:                                                                          #If last value of data is greater than las maxima
            rmax = np.append(max_array_positions[max(MAX_END - nbsym + 1, 0):], len(data_emd) - 1)[::-1]
            rmin = min_array_positions[max(MIN_END - nbsym, 0):][::-1]
            rsym = len(data_emd) - 1
    else:                                                       #if last extrema is a maxima
        if data_emd[-1] > data_emd[min_array_positions[-1]]:                          #if last data value is greater than last minima
            rmax = max_array_positions[max(MAX_END-nbsym-1,0):-1][::-1]
            rmin = min_array_positions[max(MIN_END-nbsym,0):][::-1]
            rsym = max_array_positions[-1]
        else:                                                                         #if last value of data is less than last minima
            rmax = max_array_positions[max(MAX_END-nbsym,0):][::-1]
            rmin = np.append(min_array_positions[max(MIN_END-nbsym+1,0):], len(data_emd)-1)[::-1]
            rsym = len(data_emd)-1


    #MIRROR
    tlmin = 2*x[lsym]-x[lmin]
    tlmax = 2*x[lsym]-x[lmax]
    trmin = 2*x[rsym]-x[rmin]
    trmax = 2*x[rsym]-x[rmax]



    # If mirrored points are not outside passed time range.
    if tlmin[0] > x[0] or tlmax[0] > x[0]:
        if lsym == max_array_positions[0]:
            lmax = max_array_positions[0:min(MAX_END, nbsym)][::-1]
        else:
            lmin = min_array_positions[0:min(MIN_END, nbsym)][::-1]

        lsym = 0
        tlmin = 2*x[lsym] - x[lmin]
        tlmax = 2*x[lsym] - x[lmax]

    if trmin[-1] < x[-1] or trmax[-1] < x[-1]:
        if rsym == max_array_positions[-1]:
            rmax = max_array_positions[max(MAX_END - nbsym, 0):][::-1]
        else:
            rmin = min_array_positions[max(MIN_END - nbsym, 0):][::-1]


        rsym = len(data_emd) - 1
        trmin = 2*x[rsym] - x[rmin]
        trmax = 2*x[rsym] - x[rmax]

    zlmax = data_emd[lmax]
    zlmin = data_emd[lmin]
    zrmax = data_emd[rmax]
    zrmin = data_emd[rmin]

    tmin = np.append(tlmin, np.append(x[min_array_positions], trmin))
    tmax = np.append(tlmax, np.append(x[max_array_positions], trmax))
    zmin = np.append(zlmin, np.append(data_emd[min_array_positions], zrmin))
    zmax = np.append(zlmax, np.append(data_emd[max_array_positions], zrmax))

    return tmin,tmax,zmin,zmax

fs = 500
f = 5
sample = 500
x = np.arange(sample) #x is array number of samples from 0 to maxnumber of samples

s = np.sin(2*np.pi*f*x/fs) + np.sin(10*np.pi*f*x/fs)

j = 0

component = np.zeros(shape=(x.size,10))

data_to_emd = s



while argrelmax(data_to_emd)[0].size > 3 and max(abs(data_to_emd)) > 1e-10:


    ###MAXIMA AND MINIMA CALCULATIONS
    maxima_arr = argrelmax(data_to_emd)
    minima_arr = argrelmin(data_to_emd)

    x_value_min, x_value_max, value_min, value_max = find_all_extrema(data_to_emd, x, maxima_arr[0], minima_arr[0])

    max_interpolated = interpolate.interp1d(x_value_max, value_max, kind='cubic')
    min_interpolated = interpolate.interp1d(x_value_min, value_min, kind='cubic')


    ##MEAN

    mean = (max_interpolated(x) + min_interpolated(x)) / 2

    # prints mean as before but also mean
    plt.figure()
    plt.plot(x, data_to_emd, 'r', x, max_interpolated(x), '--b', x, min_interpolated(x), '--g', x, mean, 'k')
    plt.show()

    ### ELEMENT H1 CALCULATION

    h1 = data_to_emd - mean

    h_old = h1



    sd = 1.0

    while sd > 0.2:

        ###CUBIC SPLINE CALCULATIONS

        maxima_arr = argrelmax(data_to_emd)
        minima_arr = argrelmin(data_to_emd)

        x_value_min, x_value_max, value_min, value_max = find_all_extrema(h_old, x, maxima_arr[0], minima_arr[0])

        tck_max = interpolate.splrep(x_value_max, value_max, s=0)

        tck_min = interpolate.splrep(x_value_min, value_min, s=0)

        max_interpolated = interpolate.splev(x, tck_max, der=0)

        min_interpolated = interpolate.splev(x, tck_min, der=0)  #evaluates spline for x points on y coordinates given by tck


        ##MEAN CALCULTIONS
        mean = (max_interpolated + min_interpolated)/2

        #prints mean as before but also mean
        #plt.figure()
        #plt.plot(x, s, 'r', x, max_interpolated, '--b', x, min_interpolated, '--g', x, mean, 'k')
        #plt.show()

        ### ELEMENT CALCULATION

        h_new = h_old - mean

        #prints mean as before but also mean
        plt.figure()
        plt.plot(x, h_old, 'r', x, max_interpolated, '--b', x, min_interpolated, '--g', x, h_new, 'k',x, mean, '--y')
        plt.show()

        sd = 0 #re-star SD on each loop

        for i in range (0,x.size,1):

            sd = sd + (pow(abs(h_old[i] - h_new[i]), 2))/(pow(h_old[1], 2))

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


ax1 = plt.subplot(411)
ax1.plot(s,'b')
ax2 = plt.subplot(412, sharex=ax1)
ax2.plot(component[:,0],'r')
ax3 = plt.subplot(413, sharex=ax1)
ax3.plot(component[:,1],'r')
ax4 = plt.subplot(414, sharex=ax1)
ax4.plot(data_to_emd,'r')
plt.show()