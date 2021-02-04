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
        Nosegments = int(T // ((Fs / F) * 2))  # in case freq is 40   #// is int divider
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


d, c = signal.butter(3, 50/250, 'low') 
b, a = signal.butter(3, 0.5/250, 'high')


FULLDATA_contents = sio.loadmat('Test_1mA_5Hz_ERP_PH.mat')  # FULLDATA_contents is an array of structures

DATA_struct = FULLDATA_contents['EEG']

Data = DATA_struct[0, 0]

Q = Data['Data']

tacs_t1 = 15000
tacs_t2 = tacs_t1 + 30000

Q = np.transpose(Q[1, tacs_t1:tacs_t2])

print(Q.size)


SUBJECT1_FREE_array = SMA_ERP(Q ,5 ,500)        #apply SMA TO TEST 1 DATA

#print(TEST1_FREE_array)

cProfile.run('signal.filtfilt(d, c, FREE_array, padlen=3*(max(len(d),len(c))-1)) ')

SUBJECT1_SMA_low = signal.filtfilt(d, c, SUBJECT1_FREE_array, padlen=3*(max(len(d),len(c))-1))   #FILTER TEST1 DATA AFTER SMA APPLIED
SUBJECT1_SMA = signal.filtfilt(b, a, SUBJECT1_SMA_low, padlen=3*(max(len(b),len(a))-1))



#print(ERP_TEST1_SMA)

SUBJECT2_10HZ_SMA_low = signal.filtfilt(d, c, SUBJECT2_10HZ_FREE_array, padlen=3*(max(len(d),len(c))-1))    #FILTER TEST2 DATA AFTER SMA APPLIED
SUBJECT2_10HZ_SMA = signal.filtfilt(b, a, SUBJECT2_10HZ_SMA_low, padlen=3*(max(len(b),len(a))-1))

#print(ERP_TEST2_SMA)


SHAM_SMA_low_1 = signal.filtfilt(d, c, SUBJECT1['SHAM'], padlen=3*(max(len(d),len(c))-1))           #FILTER RAW SHAM DATA TO COMPARE TO RESULTS
SHAM_SMA_1 = signal.filtfilt(b, a, SHAM_SMA_low_1, padlen=3*(max(len(b),len(a))-1))

#print(SHAM_SMA)

SHAM_SMA_low_2 = signal.filtfilt(d, c, SUBJECT2_10HZ['SHAM'], padlen=3*(max(len(d),len(c))-1))           #FILTER RAW SHAM DATA TO COMPARE TO RESULTS
SHAM_SMA_2 = signal.filtfilt(b, a, SHAM_SMA_low_2, padlen=3*(max(len(b),len(a))-1))



plot_t = np.arange(0,30000/500,1/500)

f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(plot_t,SHAM_SMA_1[0,:]/1000, '-k', label='Sham')
axarr[0].plot(plot_t,SUBJECT1_SMA[0,:]/1000, '--r', label='TEST1')
axarr[0].set_title('SUBJECT 1 and SHAM')
axarr[0].set_xlim([27,33])
axarr[0].set_ylim([-300,300])

axarr[1].plot(plot_t,SHAM_SMA_2[0,:]/1000, '-k', label='Sham')
axarr[1].plot(plot_t, SUBJECT2_10HZ_SMA[0,:]/1000, '--r', label='TEST2')
axarr[1].set_title('ERP TEST 2 and SHAM')
axarr[1].set_xlim([4,12])
axarr[1].set_ylim([-300,300])

axarr[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
axarr[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


plt.show()
