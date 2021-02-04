import numpy as np
from numpy import sqrt, square
from scipy.signal import argrelmax
from scipy.signal import argrelmin
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.io as sio
import cProfile

def _prepare_points_simple(T, S, max_pos, min_pos):

    # Find indexes of pass
    ind_min = min_pos
    ind_max = max_pos


    # Local variables
    nbsym = 2
    end_min, end_max = len(min_pos), len(max_pos)

    ####################################
    # Left bound - mirror nbsym points to the left
    if ind_max[0] < ind_min[0]:
        if S[0] > S[ind_min[0]]:
            lmax = ind_max[1:min(end_max, nbsym + 1)][::-1]
            lmin = ind_min[0:min(end_min, nbsym + 0)][::-1]
            lsym = ind_max[0]
        else:
            lmax = ind_max[0:min(end_max, nbsym)][::-1]
            lmin = np.append(ind_min[0:min(end_min, nbsym - 1)][::-1], 0)
            lsym = 0
    else:
        if S[0] < S[ind_max[0]]:
            lmax = ind_max[0:min(end_max, nbsym + 0)][::-1]
            lmin = ind_min[1:min(end_min, nbsym + 1)][::-1]
            lsym = ind_min[0]
        else:
            lmax = np.append(ind_max[0:min(end_max, nbsym - 1)][::-1], 0)
            lmin = ind_min[0:min(end_min, nbsym)][::-1]
            lsym = 0

    ####################################
    # Right bound - mirror nbsym points to the right
    if ind_max[-1] < ind_min[-1]:
        if S[-1] < S[ind_max[-1]]:
            rmax = ind_max[max(end_max - nbsym, 0):][::-1]
            rmin = ind_min[max(end_min - nbsym - 1, 0):-1][::-1]
            rsym = ind_min[-1]
        else:
            rmax = np.append(ind_max[max(end_max - nbsym + 1, 0):], len(S) - 1)[::-1]
            rmin = ind_min[max(end_min - nbsym, 0):][::-1]
            rsym = len(S) - 1
    else:
        if S[-1] > S[ind_min[-1]]:
            rmax = ind_max[max(end_max - nbsym - 1, 0):-1][::-1]
            rmin = ind_min[max(end_min - nbsym, 0):][::-1]
            rsym = ind_max[-1]
        else:
            rmax = ind_max[max(end_max - nbsym, 0):][::-1]
            rmin = np.append(ind_min[max(end_min - nbsym + 1, 0):], len(S) - 1)[::-1]
            rsym = len(S) - 1


    # Mirror points
    tlmin = (2*T[lsym] - T[lmin])
    tlmax = 2*T[lsym] - T[lmax]
    trmin = 2*T[rsym] - T[rmin]
    trmax = 2*T[rsym] - T[rmax]

    # If mirrored points are not outside passed time range.
    if tlmin[0] > T[0] or tlmax[0] > T[0]:
        if lsym == ind_max[0]:
            lmax = ind_max[0:min(end_max, nbsym)][::-1]
        else:
            lmin = ind_min[0:min(end_min, nbsym)][::-1]

        lsym = 0
        tlmin = 2*T[lsym] - T[lmin]
        tlmax = 2*T[lsym] - T[lmax]

    if trmin[-1] < T[-1] or trmax[-1] < T[-1]:
        if rsym == ind_max[-1]:
            rmax = ind_max[max(end_max - nbsym, 0):][::-1]
        else:
            rmin = ind_min[max(end_min - nbsym, 0):][::-1]

        rsym = len(S) - 1
        trmin = 2*T[rsym] - T[rmin]
        trmax = 2*T[rsym] - T[rmax]

    zlmax = S[lmax]
    zlmin = S[lmin]
    zrmax = S[rmax]
    zrmin = S[rmin]

    tmin = np.append(tlmin, np.append(T[ind_min], trmin))
    tmax = np.append(tlmax, np.append(T[ind_max], trmax))
    zmin = np.append(zlmin, np.append(S[ind_min], zrmin))
    zmax = np.append(zlmax, np.append(S[ind_max], zrmax))

    max_extrema = np.array([tmax, zmax])
    min_extrema = np.array([tmin, zmin])

    # Make double sure, that each extrema is significant
    max_dup_idx = np.where(max_extrema[0, 1:] == max_extrema[0, :-1])
    max_extrema = np.delete(max_extrema, max_dup_idx, axis=1)
    min_dup_idx = np.where(min_extrema[0, 1:] == min_extrema[0, :-1])
    min_extrema = np.delete(min_extrema, min_dup_idx, axis=1)

    return max_extrema, min_extrema

def EMD_VER1(s,MAX_ITERATIONS_NUMBER,  MAX_NUMBER_IMFS = np.inf):
    '''
    Function EMD
    Inputs: S = input signal to perform EMD on
    Outputs: plot of original signal, IMFs, residual and reconstructed signal
              RMS(s- reconstructed signal)
    '''

    s = np.reshape(s, (s.size))   #Make S a 1D array

    x = np.arange(0, s.size, 1)   #Make x axis of signal simply 1,2,3.......size of s

    j = 0 #IMFs component counter

    MAX_ITERATION = MAX_ITERATIONS_NUMBER

    component = []  #define maximum a maximum of 15 components

    data_to_emd = s

    ###ACTUAL EMD CALCULATION STARTS###
    while argrelmax(data_to_emd)[0].size > 2 and max(abs(data_to_emd)) > 1e-10:     #EMD stopping criteria based on number of extrema and

        ###DATA MAXIMA AND MINIMA CALCULATIONS
        maxima_arr = argrelmax(data_to_emd)
        minima_arr = argrelmin(data_to_emd)

        max_values, min_values = _prepare_points_simple(x, data_to_emd, maxima_arr[0], minima_arr[0])  #add values to extrapolate cubic spline

        max_interpolated = interpolate.interp1d(max_values[0], max_values[1], kind='cubic')
        min_interpolated = interpolate.interp1d(min_values[0], min_values[1], kind='cubic')

        ##MEAN
        mean = (max_interpolated(x) + min_interpolated(x)) / 2                #M(t) = (U(t)+L(t))/2

        #prints: original data, U(t), L(t), M(t)
        #plt.figure()
        #plt.plot(x, data_to_emd, 'r', label = 'Data to EMD s(t)')
        #plt.plot(x, max_interpolated(x), '--y', label = ' U(t)')
        #plt.plot(x, min_interpolated(x), '--g',  label = ' L(t)')
        #plt.plot(x, mean, 'k',  label = 'Mean')
        #plt.legend()
        #plt.show()

        ### ELEMENT H1 CALCULATION

        h1 = data_to_emd[:] - mean[:]



        #sd = np.sum((square(data_to_emd - h1)) / (square(data_to_emd)))

        sd_top = (square(data_to_emd[:]- h1[:]))
        sd_bot = (square(h1[:]))

        sd =sd_top / sd_bot

        h_old = h1

        iteration_counter = 0

        #h(k) = h(k-1)  - M(k-1) LOOP
        while np.any(sd > 0.3) and iteration_counter < MAX_ITERATION:

            ###CUBIC SPLINE CALCULATIONS
            maxima_arr_h_old = argrelmax(h_old)
            minima_arr_h_hold = argrelmin(h_old)

            max_values_h_old, min_values_h_old = _prepare_points_simple(x, h_old, maxima_arr_h_old[0], minima_arr_h_hold[0])

            max_interpolated_h_old = interpolate.interp1d(max_values_h_old[0], max_values_h_old[1], kind = 'cubic')
            min_interpolated_h_old = interpolate.interp1d(min_values_h_old[0], min_values_h_old[1], kind = 'cubic')


            ##MEAN

            mean_h_old = (max_interpolated_h_old(x) + min_interpolated_h_old(x)) / 2 #h_old M(t)

            ### NEW ELEMENT CALCULATION (h(k) = h(k-1)  - M(k-1))

            h_new = h_old - mean_h_old

            #plt.figure()
            #plt.plot(x, h_old, 'r', label = 'h_old')
            #plt.plot(x, max_interpolated_h_old(x), '--y', label = 'h_old upper')
            #plt.plot(x, min_interpolated_h_old(x), '--g',label = 'h_old lower')
            #plt.plot(x, h_new, '-b', label = 'h_new')
            #plt.plot(x, mean_h_old, 'k', label = 'h_old mean')
            #plt.legend()
            #plt.show()

            sd_top = square(h_old - h_new)
            sd_bot = square(h_new)

            sd = sd_top /sd_bot

            #sd = np.sum((square(h_old - h_new)) / (square(h_old)))   #sifiting stopping criteria

            #print(sd)

            iteration_counter = iteration_counter + 1                #sifiting counter

            h_old = h_new

        #END of IMF inner loop

        component.append(h_old)
        data_to_emd = data_to_emd - component[j]  # new data is LAST DATA - LAST EXTRACTED IMF

        if j == MAX_NUMBER_IMFS:
            return component

        print(iteration_counter)

        #PLOT of new data to emd and extracted IMF
        #lin1, axa = plt.subplots(2, sharex=True)
        #axa[0].plot(data_to_emd, 'b')
        #axa[0].set_title('Data EMD for next iteration')
        #axa[1].plot(component[j], 'r')
        #axa[1].set_title('Extracted IMF')
        #plt.show()

        j = j + 1  # look into next row of component array

    # END of EMD main loop

    residual = data_to_emd[:]

    sum_com = 0

    for m in range(j):
        sum_com = sum_com + component[m]

    initial_reconstructed = sum_com + residual

    comparison = s - initial_reconstructed

    rms = sqrt(np.mean(square(comparison)))

    f, axarr = plt.subplots(j + 3, sharex=True)
    axarr[0].plot(s, 'b')

    for k in range(1, j + 1, 1):
        axarr[k].plot(component[k - 1], 'r')

    axarr[j + 1].plot(residual, 'g')
    axarr[j + 2].plot(initial_reconstructed, 'b')
    plt.show()

    print(rms)

    if len(component) == 0:
        component.append(np.zeros(shape=(x.size)))

    return component

#Fs = 500
#f = 5
#sample = 500
#x = np.arange(sample)

#s = np.sin(2* np.pi * f * x / Fs) + np.cos(10 * np.pi * f * x / Fs) #+ np.sin(6 * np.pi * f * x / Fs)

#modes = EMD_VER1(s,20)

FULLDATA_contents = sio.loadmat('ERP_sham_test1_test2_cc.mat')  # FULLDATA_contents is an array of structures
                                                                #contains structures ERPSHAM ERPTEST1 ERPTEST2

# ALOCATE CORRECT STRUCTURES WITHIN DATA
TEST1_struct = FULLDATA_contents['ERPTEST1']    #Raw EEG+TACS DATA

# CONVERT TO NPARRAY VARIABLES
TEST1 = TEST1_struct[0, 0]

s = TEST1['E1']

s = s[0, 20000:20500]

EMD_VER1(s, 50)
#cProfile.run('EMD_VER1(s,100)')







