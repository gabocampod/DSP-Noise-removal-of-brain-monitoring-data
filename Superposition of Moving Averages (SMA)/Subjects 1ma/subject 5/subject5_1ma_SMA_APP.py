import numpy as np
import scipy.io as sio
from scipy import signal
import matplotlib.pyplot as plt


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


FULLDATA_contents = sio.loadmat('SUBJECT5_1MA_TO_PYTHON.mat')  # FULLDATA_contents is an array of structures
                                                                #contains 1 structure SUBJECT4_1MA


# ALOCATE CORRECT STRUCTURES WITHIN DATA
SUBJECT5_struct = FULLDATA_contents['SUBJECT5_1MA']  #contains EEG+Tacs values and also sham




# CONVERT TO NPARRAY VARIABLES
SUBJECT5 = SUBJECT5_struct[0, 0]


# EXAMPLES OF CALLS TO DATA STRUCTURES
# SUBJECT4['DATA_SUBJECT_4']
# SUBJECT4['SHAM_SUBJECT_4']



SUBJECT5_FREE_array = SMA_ERP(SUBJECT5['DATA_SUBJECT_5'],5,500)           #apply SMA TO  DATA
SUBJECT5_SHAM_FREE_ARRAY = SMA_ERP(SUBJECT5['SHAM_SUBJECT_5'],5,500)      #apply SMA TO  SHAM DATA



#print(TEST1_FREE_array)


SUBJECT5_SMA_low = signal.filtfilt(d, c, SUBJECT5_FREE_array, padlen=3*(max(len(d),len(c))-1))   #FILTER TEST1 DATA AFTER SMA APPLIED
SUBJECT5_SMA = signal.filtfilt(b, a, SUBJECT5_SMA_low, padlen=3*(max(len(b),len(a))-1))

#print(ERP_TEST1_SMA)


SHAM_SMA_low = signal.filtfilt(d, c, SUBJECT5_SHAM_FREE_ARRAY, padlen=3*(max(len(d),len(c))-1))           #FILTER RAW SHAM DATA TO COMPARE TO RESULTS
SHAM_SMA = signal.filtfilt(b, a, SHAM_SMA_low, padlen=3*(max(len(b),len(a))-1))

#print(SHAM_SMA)

plot_t = np.arange(0,10000/500,1/500)

fig = plt.figure();
ax = fig.add_subplot(1, 1, 1)
ax.plot(plot_t,SHAM_SMA[0,:]/1000, '-k', label='Sham')
ax.plot(plot_t,SUBJECT5_SMA[0,:]/1000, '-r', label='Subject 5 Data')
ax.set_title('Subject 5 artifact free EEG and SHAM')
ax.set_xlim([12, 17])


ax.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)

plt.show()


#print(TEST1_FREE_array)
#print(ERP_TEST1_SMA)

#working until here


sio.savemat('SUBJECT5_DATA_PYTHON.mat', {'SUBJECT5_DATA_PYTHON': SUBJECT5_SMA})
sio.savemat('SUBJECT5_SHAM_PYTHON.mat', {'SUBJECT5_SHAM_PYTHON': SHAM_SMA})
