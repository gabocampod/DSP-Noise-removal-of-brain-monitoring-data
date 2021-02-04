import numpy as np
from numpy import sqrt, square, random
from scipy.signal import argrelmax
from scipy.signal import argrelmin
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.io as sio
#from numba import jit
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

def EMD_VER1(s,MAX_ITERATIONS_NUMBER, MAX_NUMBER_IMFS = np.inf):
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

        max_interpolated = interpolate.interp1d(max_values[0], max_values[1], "cubic", fill_value='extrapolate')

        min_interpolated = interpolate.interp1d(min_values[0], min_values[1], "cubic" , fill_value='extrapolate')


        ##MEAN
        mean = (max_interpolated(x) + min_interpolated(x)) / 2      #M(t) = (U(t)+L(t))/2

        ### ELEMENT H1 CALCULATION

        h1 = data_to_emd[:] - mean[:]



        #sd = np.sum((square(data_to_emd - h1)) / (square(data_to_emd)))

        sd_top = (square(data_to_emd- h1))
        sd_bot = (square(h1))

        sd =sd_top / sd_bot

        h_old = h1

        iteration_counter = 0

        #h(k) = h(k-1)  - M(k-1) LOOP
        while np.any(sd > 0.3) and iteration_counter < MAX_ITERATION:

            ###CUBIC SPLINE CALCULATIONS
            maxima_arr_h_old = argrelmax(h_old)
            minima_arr_h_hold = argrelmin(h_old)

            max_values_h_old, min_values_h_old = _prepare_points_simple(x, h_old, maxima_arr_h_old[0], minima_arr_h_hold[0])

            max_interpolated_h_old = interpolate.interp1d(max_values_h_old[0], max_values_h_old[1], "cubic",
                                                          fill_value='extrapolate')
            min_interpolated_h_old = interpolate.interp1d(min_values_h_old[0], min_values_h_old[1], "cubic",
                                                          fill_value='extrapolate')

            ##MEAN

            mean_h_old = (max_interpolated_h_old(x) + min_interpolated_h_old(x)) / 2 #h_old M(t)

            ### NEW ELEMENT CALCULATION (h(k) = h(k-1)  - M(k-1))

            h_new = h_old - mean_h_old

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

        j = j + 1  # look into next row of component array

        if j == MAX_NUMBER_IMFS:
            component.append(data_to_emd) #return IMFS and residue
            return component

    # END of EMD main loop

    residual = data_to_emd[:]

    sum_com = 0

    for m in range(j):
        sum_com = sum_com + component[m]

    initial_reconstructed = sum_com + residual

    comparison = s - initial_reconstructed

    rms = sqrt(np.mean(square(comparison)))

    if len(component) == 0:
        component.append(np.zeros(shape=(x.size)))

    return component

def CEEMDAN_VER1(s,MAX_ITERATION_NUMBER, NR, Nstd):
    stand_deviation_s = np.std(s, ddof=1)

    s = s / stand_deviation_s

    aux = np.zeros(shape=(1, s.size))
    #iter = np.zeros(shape=(NR, int(np.round(np.log(s.size+5)))))

    white_noise = []
    modes_white_noise = []

    for i in range(NR):
        white_noise.append(random.randn(s.size))

    for i in range(NR):
        modes_white_noise.append(EMD_VER1(white_noise[i], MAX_ITERATION_NUMBER))

    #bug_modes_white_noise = np.zeros(shape=(2, s.size)) # bug_modes_white_noise[0,g] = modes_white_noise[i][0][g]
    #bug_modes_white_noise_2 = np.zeros(shape=(2, s.size)) # bug_modes_white_noise_2[0,g] = modes_white_noise_NSTD[i][0][g]
                                                            # print(bug_modes_white_noise_2[0]/bug_modes_white_noise[0])
    modes_white_noise_NSTD = modes_white_noise

    modes_white_noise_DIVI = modes_white_noise
    # []--> element on list []--> array on that element []---> element on that array
    for i in range(NR):
        for g in range(len(modes_white_noise[i][0])): # multiply all values of first IMF by Nstd
            modes_white_noise_NSTD[i][0][g] = modes_white_noise[i][0][g] * Nstd

        white_noise_stand = np.std(modes_white_noise[i][0], ddof=1)

        modes_white_noise_DIVI[i][0] = modes_white_noise_NSTD[i][0] / white_noise_stand  # divide firs IMF by first IMF's std

        si = s + modes_white_noise_DIVI[i][0]
        temp = EMD_VER1(si, MAX_ITERATION_NUMBER, 1)
        temp = temp[0]
        aux = aux + ((si - temp) / NR)

    modes = list(s - aux)
    medias = aux
    k = 0
    aux = np.zeros(shape=s.size)
    es_imf = len((EMD_VER1(medias[:], MAX_ITERATION_NUMBER, 1)))

    while es_imf > 1:  # es_imf will be one when medias has no extrema
        for i in range(NR):
            tamanio = len(modes_white_noise[i])
            if tamanio > k + 1:
                noise = modes_white_noise[i][k + 1]  # all of array number k+1 inside element 1 of list modes)
                noise = noise / np.std(noise, ddof=1)
                noise = Nstd * noise
                temp = EMD_VER1(medias[-1] + (np.std(medias[-1], ddof=1) * noise), MAX_ITERATION_NUMBER, 1)
                temp = temp[-1][:]
            else:
                temp = EMD_VER1(medias[-1], MAX_ITERATION_NUMBER, 1)  # HAPPENS LAST ITERATION
                temp = temp[-1][:]
            aux = aux + temp / NR
        modes.append(medias[-1] - aux[:])
        if type(medias) == np.ndarray:
            medias = list(medias)
        medias.append(aux)
        aux = np.zeros(shape=s.size)
        k = k + 1
        es_imf = len(EMD_VER1(medias[-1], MAX_ITERATION_NUMBER, 1))
    #end CEEMDAN calculation
    modes.append(medias[-1])
    modes_final = [[D * stand_deviation_s for D in y] for y in modes]  # modes = modes*desvio s modes is 2d array

    sum_com = np.zeros(shape=(1, s.size))

    for m in range(len(modes_final)):
        sum_com[0] += modes_final[m][:]
    #'''''
    f, axam = plt.subplots(len(modes_final) + 2, sharex=True)
    axam[0].plot(s, 'b')

    for p in range(1, len(modes_final) + 1, 1):
        axam[p].plot(modes_final[p - 1], 'r')

    axam[len(modes_final) + 1].plot(sum_com[0])
    plt.show()
    #'''''
    return modes_final

####MAIN LOOP####

FULLDATA_contents = sio.loadmat('ERP_sham_test1_test2_cc.mat')  # FULLDATA_contents is an array of structures
                                                                #contains structures ERPSHAM ERPTEST1 ERPTEST2

# ALOCATE CORRECT STRUCTURES WITHIN DATA
TEST1_struct = FULLDATA_contents['ERPTEST1']    #Raw EEG+TACS DATA

# CONVERT TO NPARRAY VARIABLES
TEST1 = TEST1_struct[0, 0]

s = TEST1['E1']


s = s[0,5000:25000]


#final_modes = CEEMDAN_VER1(s,5, 50, 0.2)
cProfile.run('CEEMDAN_VER1(s,5,100,0.02)')

#sio.savemat('CEEMDAN_2000_50_500_02_RESULTS.mat', {'CEEMDAN_2000_50_500_02': final_modes})
