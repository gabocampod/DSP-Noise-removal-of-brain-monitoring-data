import numpy as np

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
            rsym = min_array_positions[-1]
        else:                                                                          #If last value of data is greater than last maxima
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