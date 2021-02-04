import numpy as np


def _prepare_points_simple(T, S, max_pos, min_pos):

    # Find indexes of pass
    ind_min = np.array([np.nonzero(T == t)[0] for t in min_pos]).flatten()
    ind_max = np.array([np.nonzero(T == t)[0] for t in max_pos]).flatten()

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

    # In case any array missing
    if not lmin.size: lmin = ind_min
    if not rmin.size: rmin = ind_min
    if not lmax.size: lmax = ind_max
    if not rmax.size: rmax = ind_max

    # Mirror points
    tlmin = 2*T[lsym] - T[lmin]
    tlmax = 2*T[lsym] - T[lmax]
    trmin = 2*T[rsym] - T[rmin]
    trmax = 2*T[rsym] - T[rmax]

    # If mirrored points are not outside passed time range.
    if tlmin[0] > T[0] or tlmax[0] > T[0]:
        if lsym == ind_max[0]:
            lmax = ind_max[0:min(end_max, nbsym)][::-1]
        else:
            lmin = ind_min[0:min(end_min, nbsym)][::-1]

        if lsym == 0:
            raise Exception('Left edge BUG')

        lsym = 0
        tlmin = 2*T[lsym] - T[lmin]
        tlmax = 2*T[lsym] - T[lmax]

    if trmin[-1] < T[-1] or trmax[-1] < T[-1]:
        if rsym == ind_max[-1]:
            rmax = ind_max[max(end_max - nbsym, 0):][::-1]
        else:
            rmin = ind_min[max(end_min - nbsym, 0):][::-1]

        if rsym == len(S) - 1:
            raise Exception('Right edge BUG')

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

    # Make double sure, that each extremum is significant
    max_dup_idx = np.where(max_extrema[0, 1:] == max_extrema[0, :-1])
    max_extrema = np.delete(max_extrema, max_dup_idx, axis=1)
    min_dup_idx = np.where(min_extrema[0, 1:] == min_extrema[0, :-1])
    min_extrema = np.delete(min_extrema, min_dup_idx, axis=1)

    return max_extrema, min_extrema