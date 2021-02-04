import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

def SMA_REAL_TIME(INPUT,currentfreq,SamplingFreq):

    data = np.zeros(shape=(1, INPUT.size))  # i.e #data is horizontal [0 0 0 0 ]

    data = INPUT

    datareshape = np.reshape(data, (data.size, 1))  # now data is vertical [0 0 0 ]´

    T = len(datareshape)  # lenght of DATA
    F = currentfreq       # Current frequency (could be asked to user)
    Fs = SamplingFreq     # sampling frequency (could be asked to user)

    if F == 40:
        Nosegments = (T // ((Fs / F) * 2))  # in case freq is 40   #// is int divider
    else:
        Nosegments = (T // ((Fs // F) * 1))  # for 5h and 10hz

    AdjSegments = 10
    SegmentCount = AdjSegments + 1
    timestep = (T // Nosegments)
    Start = int(0.4 * Fs)
    StartSeg = int(((T - Start) // timestep) + 1)

    # Assuming 1 channel

    # used variables
    segments = np.zeros(shape=(timestep, Nosegments))  # i.e timestep x Nosegments array
    CurrentSegments = np.zeros(shape=(timestep, SegmentCount))
    ScaledArtefact = np.zeros(shape=(timestep, (Nosegments - StartSeg) + 1))
    AvgArtefact = np.zeros(shape=(Start, 1))
    Signal = np.zeros(shape=(1, Start))
    # start of code

    t1 = 0
    t2 = timestep

    for j in range(Nosegments):
        segments[:, j] = datareshape[t1:t2, 0]
        t1 = t1 + timestep
        t2 = t2 + timestep

    ##print(segments)
    a = 0
    for k in range(StartSeg, Nosegments + 1, 1):
        S1 = StartSeg - AdjSegments
        S2 = StartSeg
        CurrentSegments[:, 0:SegmentCount] = segments[:, S1 - 1:S2]
        ScaledArtefact[:, a] = CurrentSegments.mean(axis=1)
        a = a + 1

    AvgArtefact = np.reshape(ScaledArtefact, (Start, 1), order='F')
    data = np.reshape(data, (1, T))
    AvgArtefact = np.reshape(AvgArtefact, (1, Start))
    Signal = data[0, (T - Start):T] - AvgArtefact

    #print(Signal)

    return Signal

    '''
                            FOR DEBUGGING 
    print(Nosegments)
    print(timestep)
    print(Start)
    print(StartSeg)
    print(segments.shape)
    print(CurrentSegments.shape)  ##ok
    print(ScaledArtefact.shape)
    print(AvgArtefact.shape)
    print(Signal.shape)
    '''




FULLDATA_contents = sio.loadmat('DATEtime.mat')  # FULLDATA_contents is an array of structures
                                                 #DATAtime contains four structures ALPHA1 ERP1  SHAMALPHA1 AND SHAMERP1

# ALOCATE CORRECT STRUCTURES WITHIN DATA
ALPHA_struct = FULLDATA_contents['ALPHA1']  # ALPHA_STRUCTS is an array of structures with 1 element
ERP_struct = FULLDATA_contents['ERP1']
SHAMALPHA_struct = FULLDATA_contents['SHAMALPHA1']
SHAMERP_struct = FULLDATA_contents['SHAMERP1']

                                            #Each of this is a 1x1 struct with two fields A1 and TIME

# CONVERT TO NPARRAY VARIABLES
ALPHA_VALandTIME = ALPHA_struct[0, 0]  # ALPHA_VALandTIME is an array of two elements
ERP_VALandTIME = ERP_struct[0, 0]
SHAMALPHA_VALandTIME = SHAMALPHA_struct[0, 0]
SHAMERP_VALandTIME = SHAMERP_struct[0, 0]

                                        #now we are saying alphaValandTIME is element 1 of alpha struct in other words we copy all of
                                        #one into the other but the new one is easier to call

# EXAMPLES OF CALLS TO DATA STRUCTURES
# ALPHA_VALandTIME['A1']
# ALPHA_VALandTIME['TIME']

##print(ALPHA_VALandTIME['A1'].shape)
##print(ALPHA_VALandTIME['A1'])
##print(type(ALPHA_VALandTIME['A1']))

ALPHA_FREE_array = SMA_REAL_TIME(ALPHA_VALandTIME['A1'],5,500)

ERP_FREE_array = SMA_REAL_TIME(ERP_VALandTIME['E1'],5,500)

datatryyy = np.zeros(shape=(1,500))

print(datatryyy.shape)

fig, axs = plt.subplots(2, 1, constrained_layout=True)
axs[0].stem(ALPHA_FREE_array[0]/1000)
axs[0].set_title('ALPHA EEG')
axs[0].set_xlabel('Time')
axs[0].set_ylabel('EEG data')
fig.suptitle('ALPHA and ERP results', fontsize=16)

axs[1].stem(ERP_FREE_array[0]/1000)
axs[1].set_title('ERP EEG')
axs[1].set_xlabel('Time')
axs[1].set_ylabel('EEG data')


plt.show()
