import numpy as np
import scipy.io as sio
from scipy import signal
import matplotlib.pyplot as plt
import cProfile


def SMA_ERP(INPUT,currentfreq,SamplingFreq):

    data = np.zeros(shape=(1, INPUT.size))  # i.e #data is horizontal [0 0 0 0 ]

    data = INPUT

    datareshape = np.reshape(data, (data.size, 1))  # now data is vertical [0 0 0 ]´

    T = len(datareshape)  # length of DATA
    if np.mod(T,2)==1:
        T = T-1
    F = currentfreq  # Current frequency (could be asked to user)
    Fs = SamplingFreq  # sampling frequency (could be asked to user)

    TIME = T / Fs
    my_new_t = np.arange(1 / Fs, TIME, 1 / Fs)

    if F == 40:
        Nosegments = (T // ((Fs / F) * 2))  # in case freq is 40   #// is int divider
        AdjSegments = np.floor(0.05 * Nosegments)  # ADJ is an integer still
    else:
        Nosegments = (T // ((Fs // F) * 1))  # for 5h and 10hz
        AdjSegments = np.floor(0.05 * Nosegments)

    if np.mod(AdjSegments, 2) == 1:  # if adj segments is odd then:
        AdjSegments = AdjSegments - 1;  # make adjsegments even

    SegmentCount = int(AdjSegments + 1)  # still an int
    timestep = (T // Nosegments)
    # Assuming 1 channel

    # used variables
    segments = np.zeros(shape=(timestep, Nosegments))  # i.e timestep x Nosegments array
    CurrentSegments = np.zeros(shape=(timestep, SegmentCount))  # should still be true
    ScaledArtefact = np.zeros(shape=(timestep, Nosegments))
    AvgArtefact = np.zeros(shape=(T, 1))
    Signal = np.zeros(shape=(1, T))
    # start of code

    t1 = 0
    t2 = timestep

    for j in range(Nosegments):  # Has not changed
        segments[:, j] = datareshape[t1:t2, 0]
        t1 = t1 + timestep
        t2 = t2 + timestep

    for k in range(1, Nosegments + 1, 1):  # I added that + 1
        if k <= (AdjSegments / 2):  # remember ADJ segments is even so this always gives int
            S1 = 0
            S2 = AdjSegments
            CurrentSegments[:, 0:SegmentCount] = segments[:, 0:SegmentCount]
        else:
            if k >= (Nosegments - (AdjSegments / 2)):
                CurrentSegments[:, 0:SegmentCount] = segments[:, (Nosegments - SegmentCount):Nosegments]
            else:
                CurrentSegments[:, 0:SegmentCount] = segments[:,(int(k - AdjSegments / 2) - 1):int((k + AdjSegments / 2))]

        ScaledArtefact[:, k - 1] = CurrentSegments.mean(axis=1)

    AvgArtefact = np.reshape(ScaledArtefact, (T, 1), order='F')  # DEFINE AVGARTEFACT

    AvgArtefact = np.reshape(AvgArtefact, (1, T))  # TURN INTO HORIZONTAL VECTOR

    Signal = data[0,0:T] - AvgArtefact

    # print(Signal)

    return Signal

    '''
    print(Nosegments) 
    print(timestep)
    print(AdjSegments)
    print(SegmentCount)
    print(segments.shape)
    print(CurrentSegments.shape)  ##ok
    print(ScaledArtefact.shape)
    print(AvgArtefact.shape)
    print(Signal.shape)
    '''

d, c = signal.butter(3, 50/250, 'low')
b, a = signal.butter(3, 0.5/250, 'high')


FULLDATA_contents = sio.loadmat('ERP_sham_test1_test2_cc.mat')  # FULLDATA_contents is an array of structures
                                                            #contains four structures ERPSHAM ERPTEST1 ERPTEST2


# ALOCATE CORRECT STRUCTURES WITHIN DATA
SHAM_ERP_struct = FULLDATA_contents['ERPSHAM']  #Raw Sham Data
TEST1_struct = FULLDATA_contents['ERPTEST1']    #Raw EEG+TACS DATA
TEST2_struct = FULLDATA_contents['ERPTEST2']    #Raw EEG+TACS DATA
                                                #Each of this is a 1x1 struct with two fields


# CONVERT TO NPARRAY VARIABLES
SHAM = SHAM_ERP_struct[0, 0]
TEST1 = TEST1_struct[0, 0]
TEST2 = TEST2_struct[0, 0]
                                                # now we are saying TEST1 is element 1 of alpha struct in other words we copy all of
                                                # one into the other but the new one is easier to call

# EXAMPLES OF CALLS TO DATA STRUCTURES
# TEST1['E1']
# ALPHA_VALandTIME['FREE_EEG_1']



TEST1_FREE_array = SMA_ERP(TEST1['E1'],5,500)           #apply SMA TO TEST 1 DATA
TEST2_FREE_array = SMA_ERP(TEST2['E2'],5,500)           #apply SMA

#cProfile.run('SMA_ERP(TEST2["E2"],5,500)')

#print(TEST1_FREE_array)


ERP_TEST1_SMA_low = signal.filtfilt(d, c, TEST1_FREE_array, padlen=3*(max(len(d),len(c))-1))   #FILTER TEST1 DATA AFTER SMA APPLIED
ERP_TEST1_SMA = signal.filtfilt(b, a, ERP_TEST1_SMA_low, padlen=3*(max(len(b),len(a))-1))

#print(ERP_TEST1_SMA)



ERP_TEST2_SMA_low = signal.filtfilt(d, c, TEST2_FREE_array, padlen=3*(max(len(d),len(c))-1))    #FILTER TEST2 DATA AFTER SMA APPLIED
ERP_TEST2_SMA = signal.filtfilt(b, a, ERP_TEST2_SMA_low, padlen=3*(max(len(b),len(a))-1))

#print(ERP_TEST2_SMA)



SHAM_SMA_low = signal.filtfilt(d, c, SHAM['ES'], padlen=3*(max(len(d),len(c))-1))           #FILTER RAW SHAM DATA TO COMPARE TO RESULTS
SHAM_SMA = signal.filtfilt(b, a, SHAM_SMA_low, padlen=3*(max(len(b),len(a))-1))

#print(SHAM_SMA)



plot_t = np.arange(0,30000/500,1/500)

f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(plot_t,SHAM_SMA[0,:]/1000, '-k', label='Sham')
axarr[0].plot(plot_t,ERP_TEST1_SMA[0,:]/1000, '--r', label='TEST1')
axarr[0].set_title('ERP TEST 1 and SHAM')
axarr[0].set_xlim([27,33])
axarr[0].set_ylim([-100,100])

axarr[1].plot(plot_t,SHAM_SMA[0,:]/1000, '-k', label='Sham')
axarr[1].plot(plot_t, ERP_TEST2_SMA[0,:]/1000, '--r', label='TEST2')
axarr[1].set_title('ERP TEST 2 and SHAM')
axarr[1].set_xlim([27,33])
axarr[1].set_ylim([-100,100])
plt.xlabel('time (s)')

axarr[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
axarr[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


plt.show()


#print(TEST1_FREE_array)
#print(ERP_TEST1_SMA)

#working until here

sio.savemat('ERP_TEST1_BEFORE_FILTERS.mat', {'TEST1_BEFORE': TEST1_FREE_array})
sio.savemat('ERP_TEST1_AFTER_FILTERS.mat', {'TEST1_AFTER': ERP_TEST1_SMA})
sio.savemat('ERP_TEST2_BEFORE_FILTERS.mat', {'TEST2_BEFORE': TEST2_FREE_array})
sio.savemat('ERP_TEST2_AFTER_FILTERS.mat', {'TEST2_AFTER': ERP_TEST2_SMA})

sio.savemat('ERP_PYTHON_SHAM_AFTER_FILTER.mat', {'ERP_SHAM_AFTER': SHAM_SMA})